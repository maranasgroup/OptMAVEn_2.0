__doc__ = """ This module provides a consistent interface to the MAPs
database. """

import cPickle as pkl
import itertools
import os
import re
import sys

import numpy as np

import standards
sys.path.append(standards.CplexDirectory)
import cplex

coords_files = dict()
for gap in standards.MapsGaps:
    coords_files[gap] = os.path.join(standards.MapsDirectory, "coords_{}.csv".format(
            gap))

sep = "_"
light_chains = standards.MapsLightChains
heavy_chains = standards.MapsHeavyChains
chains = standards.MapsChains
regions = standards.MapsRegions
cdrs = standards.MapsCdrs
cdr_directories = [os.path.join(standards.MapsDirectory, cdr) for cdr in cdrs]
parts = {os.path.splitext(part)[0]: os.path.join(cdr, part) for cdr in cdr_directories for part in os.listdir(cdr)}
part_cdr = dict()
part_number = dict()
part_pattern = re.compile("([A-Z3]+)_([0-9]+)")
for part in parts:
    match = part_pattern.match(part)
    if not match:
        raise ValueError("Bad part name: {}".format(part))
    cdr, number = match.groups()
    if cdr not in cdrs:
        raise ValueError("Bad part cdr: {}".format(cdr))
    try:
        int(number)
    except ValueError:
        raise ValueError("Bad part number: {}".format(number))
    part_cdr[part] = cdr
    part_number[part] = number


coords = dict()
for gap in standards.MapsGaps:
    coords[gap] = {line.split(",")[0]: np.array(map(float, line.split(",")[1:]))
            for line in open(coords_files[gap])}


def validate_part(part):
    if part not in parts:
        raise ValueError("Bad MAPs part: {}".format(part))


def join(cdr, number, check_new_part=True):
    """ Join a tuple X, n into a part X_n """
    part = "{}{}{}".format(cdr, sep, number)
    if check_new_part and part not in parts:
        raise ValueError("Bad MAPs part: {}".format(part))
    return part


def split(part, check_part=True):
    """ Split a part X_n into a tuple X, n """
    if check_part and part not in parts:
        raise ValueError("Bad MAPs part: {}".format(part))
    return tuple(part.split(sep))


def join_cdr(cdr, region):
    """ Join a chain and region into a CDR. """
    cdr = "{}{}".format(cdr, region)
    if cdr not in cdrs:
        raise ValueError("Bad CDR: {}".format(cdr))
    return cdr 


def split_cdr(item):
    """ Split a part or CDR into the chain and region. """
    try:
        item, number = split(item)
    except ValueError:
        pass
    if item in cdrs:
        return (item[0], item[1:])
    else:
        raise ValueError("Bad MAPs item: {}".format(item))


def get_part_cdr(part):
    return split(part)[0]


def get_part_number(part):
    return split(part)[1]


def get_chain(item):
    """ Get the chain (H, K, or L) of a part or cdr. """
    return split_cdr(item)[0]


def get_region(item):
    return split_cdr(item)[1]


def translate_chain(item, new_chain, check_new_part=True):
    """ Convert the chain of one MAPs part into another chain. """
    if new_chain not in chains:
        raise ValueError("Bad chain: {}".format(new_chain))
    try:
        item, number = split(item)
    except ValueError:
        number = None
    new_cdr = join_cdr(new_chain, get_region(item))
    if number is None:
        return new_cdr
    else:
        return join(new_cdr, number, check_new_part)


def translate_chain_namesake(item):
    """ Convert the chain of a MAPs part into its namesake value: H -> H, K -> L, L -> L. """
    chain = get_chain(item)
    if chain in standards.MapsHeavyChains:
        return translate_chain(item, standards.MapsNamesakeHeavy)
    elif chain in standards.MapsLightChains:
        return translate_chain(item, standards.MapsNamesakeLight)
    else:
        raise ValueError("Bad chain: {}".format(chain))


def get_cdr_names(light_chain):
    """ Get the six CDRs for a given light chain. """
    if light_chain not in light_chains:
        raise ValueError("Bad light chain: {}".format(light_chain))
    return [cdr for cdr in cdrs if get_chain(cdr) in heavy_chains + [light_chain]]


def get_integer_cuts(expanded_format=True):
    """ Load all of the integer cuts of MAPs parts that clash. """
    cuts = list()
    with open(standards.MapsIntegerCutsFile) as f:
        for line in f:
            # Skip blank lines.
            if line.strip() == "":
                continue
            cut = line.split()
            try:
                cdr1, num1, cdr2, num2 = cut
                part1 = join(cdr1, num1)
                part2 = join(cdr2, num2)
            except ValueError:
                raise ValueError("Bad integer cut: {}".format(line))
            if expanded_format:
                cuts.append(cut)
            else:
                cuts.append((part1, part2))
    return cuts


maps_clashes = set(get_integer_cuts(expanded_format=False))
def clashing(part1, part2, error=False):
    """ Check if two MAPs parts clash. """
    validate_part(part1)
    validate_part(part2)
    if (part1, part2) in maps_clashes or (part2, part1) in maps_clashes:
        if error:
            raise ValueError("MAPs parts clash: {}, {}".format(part1, part2))
        else:
            return True
    else:
        return False


def select_parts(energies, clash_cuts, solution_cuts):
    """Use CPLEX to select an optimal combination of MAPs parts"""
    # Make the model and get a dictionary containing the different kinds of
    # parts
    model, parts = make_optmaven_selector(energies, clash_cuts, solution_cuts)
    # Solve the model
    model.solve()
    # Get the status of the solution
    status = model.solution.get_status()
    if status != 101 and status != 102:
        raise ValueError("CPLEX exited with status {}".format(status))
    # Get the objective value
    objective = model.solution.get_objective_value()
    # Store the solution here
    solution = {}
    for part in parts:
        for i in parts[part]:
            # The name of this variable
            name = "X_" + part + "_" + str(i)
            # The value of the variable
            x = model.solution.get_values(name)
            if x > 0.01:
                if part not in solution:
                    solution[part] = i
                else:
                    raise ValueError("The CPLEX solution has multiple parts selected for {}".format(part))
    return solution, objective


def make_optmaven_selector(energies, clash_cuts, solution_cuts):
    """Make a CPLEX model to select an optimal combination of MAPs parts"""
    # Generate a CPLEX model
    model = cplex.Cplex()
    # Set it up as an MILP minimizing the energy
    model.set_problem_type(cplex.Cplex.problem_type.MILP)
    model.objective.set_sense(model.objective.sense.minimize)
    # Make a list of the parts, divided by region
    parts = {}
    for data in energies:
        if data[0] not in parts:
            parts[data[0]] = []
        parts[data[0]].append(data[1])
    # Determine the types of domains that are being created
    domains = []
    for part in parts:
        if part.startswith("H") and "H" not in domains:
            domains.append("H")
        elif part.startswith(("K", "L")) and "L" not in domains:
            domains.append("L")
    # Create the objective function
    # Store the variables here
    objV = []
    # and their coefficients (i.e. energies) here
    objC = []
    # go through the energies
    for data in energies:
        objV.append("X_" + data[0] + "_" + str(data[1]))
        objC.append(data[2])
    model.variables.add(names = objV, obj = objC, lb = [0] * len(objV), ub = \
                        [1] * len(objV), types = ['B'] * len(objV))
    # Now create the restraints. First, make sure exactly one heavy V* structure
    # is selected
    if "H" in domains:
        vars = []
        for i in parts["HV"]:
            vars.append("X_HV_" + str(i))
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, [1] * \
                                     len(vars))], senses = ["E"], rhs = [1])
        vars = []
        for i in parts["HCDR3"]:
            vars.append("X_HCDR3_" + str(i))
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, [1] * \
                                     len(vars))], senses = ["E"], rhs = [1])
        vars = []
        for i in parts["HCDR3"]:
            vars.append("X_HCDR3_" + str(i))
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, [1] * \
                                     len(vars))], senses = ["E"], rhs = [1])
        vars = []
        for i in parts["HJ"]:
            vars.append("X_HJ_" + str(i))
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, [1] * \
                                     len(vars))], senses = ["E"], rhs = [1])

    # Now make sure only one light V* is selected
    if "L" in domains:
        vars = []
        for p in ["KV", "LV"]:
            if p not in parts:
                continue
            for i in parts[p]:
                vars.append("X_" + p + "_" + str(i))
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, [1] * \
                                     len(vars))], senses = ["E"], rhs = [1])

        # Now make sure exactly one CDR3 and J* structure is selected for each type
        # of variable domain
        for L in ["K", "L"]:
            if L + "V" not in parts:
                continue
            # Store the V* info here
            Vvars = []
            Vcoefs = []
            for i in parts[L + "V"]:
                Vvars.append("X_" + L + "V_" + str(i))
                Vcoefs.append(-1)
            # Do this for each of the other two regions of structure
            for R in ['CDR3', 'J']:
                if L + R not in parts:
                    continue
                vars = []
                coefs = []
                for i in parts[L+R]:
                    vars.append("X_" + L + R + "_" + str(i))
                    coefs.append(1)
                # Include the V* information in the restraint
                vars.extend(Vvars)
                coefs.extend(Vcoefs)
                # Make the restraint
                model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, \
                                            coefs)], senses = ["E"], rhs = [0])
    # Include the integer cuts on structure combinations
    for cut in clash_cuts:
        if cut[0] in parts and cut[1] in parts[cut[0]] and \
                cut[2] in parts and cut[3] in parts[cut[2]]:
            # Say that those two cannot be selected together
            vars = ["X_" + cut[0] + "_" + str(cut[1]),
                    "X_" + cut[2] + "_" + str(cut[3])]
            coefs = [1, 1]
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, \
                                         coefs)], senses = ["L"], rhs = [1])
    # Include the integer cuts from previous optimized solutions for finding the next optimized one
    # First delete the HJ or LJ or KJ from the solution cuts because they almost contribute nothing to the energy
    # Otherwise you can get almost the same parts solution except the 'J' parts are different
    newSolutionCuts = []
    for cut in solution_cuts:
        index = 0
        newCut = []
        length = len(cut)
        while(index < length-1):
            if not cut[index][1] == 'J':
                newCut.append(cut[index])
                newCut.append(cut[index + 1])
            index += 2
        newSolutionCuts.append(newCut)
    for cut in newSolutionCuts:
        if cut[0] in parts and cut[1] in parts[cut[0]] and \
           cut[2] in parts and cut[3] in parts[cut[2]] and \
           cut[4] in parts and cut[5] in parts[cut[4]] and \
           cut[6] in parts and cut[7] in parts[cut[6]]:
            vars = ["X_" + cut[0] + "_" + str(cut[1]), \
                    "X_" + cut[2] + "_" + str(cut[3]), \
                    "X_" + cut[4] + "_" + str(cut[5]), \
                    "X_" + cut[6] + "_" + str(cut[7])]
            coefs = [1, 1, 1, 1]
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, \
                                         coefs)], senses = ["L"], rhs = [3])

    # Return the model and the dictionary of parts
    return model, parts



def get_coordinates(part, gap):
    return coords[gap][part]


def MAPs_part_clashes(cut_file, clash_cutoff=1.0):
    """ This function is only used for the pre-processing step of
    determining which pairs of MAPs parts clash sterically and
    developing integer cuts. """
    parts = get_parts()
    # Choose all combinations of categories.
    # Do not check for clashes between parts of the same category.
    finished = True
    text = ""
    for cat1, cat2 in itertools.combinations(parts, 2):
        # Also, do not check for clashes between one kappa and one lambda part.
        if "".join(sorted(cat1[0] + cat2[0])) != "KL":
            result = "{}_{}_{}.txt".format(cut_file, cat1, cat2)
            if os.path.isfile(result):
                with open(result) as f:
                    text += f.read().strip() + "\n"
            else:
                script = "{}_{}_{}.sh".format(cut_file, cat1, cat2)
                command = "python {}/Make_Integer_Cuts.py {} {} {} {}".format(standards.MapsDirectory, cat1, cat2, result, clash_cutoff)
                submitter.experiment_script(None, script, command, self_destruct=True)
                finished = False
    if finished:
        with open(cut_file + ".txt", "w") as f:
            f.write(text)
        for cat1, cat2 in itertools.combinations(parts, 2):
            if "".join(sorted(cat1[0] + cat2[0])) != "KL":
                result = "{}_{}_{}.txt".format(cut_file, cat1, cat2)
                try:
                    os.remove(result)
                except OSError:
                    pass
