""" This module defines the Experiment class of OptMAVEn. """

import argparse
from collections import defaultdict, OrderedDict
import cPickle as pkl
import csv
import os
import shutil
import sys
import tempfile
import time
import traceback

from Bio.PDB import Residue, Selection
from Bio.PDB.PDBIO import Select
from Bio import SeqIO
import numpy as np

import benchmarking
import charmm
from console import clear, disp
import klaus
import kmeans
import maps
import molecules
import standards
import submitter
import user_input


class Experiment(object):
    """ An Experiment is the class that orchestrates all tasks in OptMAVEn. """
    def __init__(self, args=None):
        # When an Experiment is created, its first purpose is to initialize itself.
        self.purpose = "Initialization"
        # Import arguments from the command line.
        if args is None:
            args = ExperimentArgParser().parse_args()
        self.args = args
        clear()
        # Let the user name the Experiment using a name that does not already exist.
        do = True
        while do:
            name = self.ask("name", valid_path=True)
            directory = os.path.join(standards.ExperimentsDirectory, name)
            if os.path.isdir(directory):
                disp("There is already an Experiment named {}.".format(name))
            else:
                do = False
        # Set up the names of the Experiment and its directories.
        self.name = name
        self.name_contig = "_".join(name.split())
        self.directory = directory  # main Experiment directory
        os.mkdir(self.directory)
        self.report_directory()
        self.temp = os.path.join(self.directory, ".temp")  # directory for temporary files: deleted at the end unless --keeptemp is given as an argument
        self.get_temp()
        self.file = os.path.join(self.temp, "{}.pickle".format(self.name_contig))  # pickle file to store the Experiment object
        self.warnings = os.path.join(self.directory, "warnings.txt")  # file of warnings
        self.errors = os.path.join(self.directory, "errors.txt")  # file of errors that may occur
        self.structure_directory = os.path.join(self.directory, "structures")
        os.mkdir(self.structure_directory)
        # Ask whether to turn on benchmarking for the Experiment.
        self.benchmarking = user_input.get_yn("Would you like to turn on benchmarking? ")
        if self.benchmarking:
            self.benchmark_directory = os.path.join(self.get_temp(), "benchmarking")
            os.mkdir(self.benchmark_directory)
            # Record the drive usage at the very beginning.
            self.add_drive_usage()

    def save(self):
        """ Save the current state of the Experiment to the pickle file. """
        with open(self.file, "w") as f:
            pkl.dump(self, f)

    def submit(self, args=None, jobs=None, options=None, queue=True):
        """ Submit a job or batch of jobs to the PBS queue. """
        if args is None:
            args = [""]
        if jobs is None:
            # If there are no jobs specified, then just run as one job (i.e. not as a job array).
            # Create a file for the script that will be submitted.
            handle, file_name = tempfile.mkstemp(prefix="{}_".format(self.name_contig), suffix=".sh", dir=self.temp)
            os.close(handle)
            # Create the command within the aforementioned script that will run whatever task needs to be submitted.
            # This command takes the form '/path/to/python /path/to/experiment.py /path/to/experiment_name.pickle [arg1] [arg2] ...'
            # When experiment.py is run directly (i.e. __name__ == '__main__') it unpickles and runs the Experiment object that follows in the command.
            command = "{} {} {} {}".format(standards.PythonCommand, os.path.realpath(__file__), self.file, " ".join(map(str, args)))
            if self.benchmarking:
                # Make a time file in which to store the time that this job took.
                time_file = self.make_time_file()
                self.add_time_file(time_file)
            else:
                time_file = None
            # Submit the job.
            submitter.submit(file_name, command, self.walltime, options=options, queue=queue, purpose=self.purpose, time_file=time_file)
        else:
            # If jobs have been specified, then use the PbsBatchSubmitter class to handle the jobs.
            s = submitter.PbsBatchSubmitter(self)
            s.submit(standards.PythonCommand, [os.path.realpath(__file__), self.file], jobs)

    def save_and_submit(self, queue=True):
        """ Save before submitting to ensure that the Experiment is run from the proper state. """
        self.save()
        self.submit(queue=queue)

    def add_benchmark(self, task):
        """ Add a pickled benchmarking Task object to the benchmark directory. These objects will be collected at the end of the Experiment. """
        handle, file_name = tempfile.mkstemp(dir=self.benchmark_directory, suffix=".pickle")
        os.close(handle)
        with open(file_name, "w") as f:
            pkl.dump(task, f)

    def add_time_file(self, _file, status_offset=0):
        """ Add a TimeFile record to the experiment. """
        task, purpose = self.get_task(self.status + status_offset)
        task = benchmarking.TimeFile(_file, purpose)
        self.add_benchmark(task)

    def add_drive_usage(self, _file=True):
        """ Add a DriveUsage record to the experiment. """
        du = benchmarking.DriveUsage(self)
        if _file:
            self.add_benchmark(du)
        return du

    def ask(self, attribute, number=False, valid_path=False):
        """ Ask for an attribute of the experiment.
        number: must the attribute be numeric?
        valid_path: must the attribute be a valid path name?
        """
        try:
            name = self.name
        except AttributeError:
            name = "this Experiment"
        prompt = "Please enter the {} of {}: ".format(attribute, name)
        do = True
        while do:
            if number:
                answer = user_input.get_number(prompt)
            else:
                answer = user_input.get(prompt)
            if not valid_path or standards.is_path_component(answer):
                do = False
            else:
                disp("The {} of {} must be a valid component of a path.".format(attribute, name))
        return answer
    
    def report_directory(self):
        """ Print the location of the Experiment results. """
        try:
            disp("The results of {} will be located in the following directory:\n{}".format(self.name, self.directory))
        except AttributeError:
            disp("This Experiment has no directory.")

    def get_temp(self):
        """ Return the directory of temporary files. If it does not exist, create it. This is to avoid returning a nonexistent directory. """
        if not os.path.isdir(self.temp):
            os.mkdir(self.temp)
        return self.temp

    def make_time_file(self):
        """ Make a temporary file in which to store a time for a TimeFile object. """
        handle, name = tempfile.mkstemp(dir=self.get_temp(), prefix="time_", suffix=".txt")
        os.close(handle)
        return name

    def document_error(self, message):
        """ Record an error that has occurred in the errors.txt file. """
        try:
            # Try to write the error in the proper file.
            with open(self.errors, "a") as f:
                f.write(str(message))
        except OSError as e:
            # If that fails for some reason, then write in a file called Optmaven_Experiment_errors.txt in the current working directory.
            error_file = "Optmaven_Experiment_errors.txt"
            with open(error_file, "a") as f:
                f.write("{}\n{}".format(message, e.message))

    def run(self, args=None):
        """ Run whatever task is indicated by the Experiment's present status. """
        try:
            try:
                task, self.purpose = self.get_task(self.status)
            except (IndexError, TypeError):
                raise ValueError("Bad Experiment status: {}".format(self.status))
            else:
                task(args)
        except Exception as e:
            # Catch all errors that occur while running a task so that they can be documented.
            # Document their tracebacks as well to be as informative as possible.
            tb = traceback.format_exc()
            self.document_error(tb)

    def run_next(self, args=None):
        """ Run the next task by incrementing the status and then running. """
        self.change_status()
        self.run(args)

    def get_molecule(self, name, prompt):
        """ Ask the user for the name of a molecule file and load the molecule. """
        do = True
        while do:
            # Ask for the file name. If it doesn't exist but is a valid PDB ID, then fetch the PDB.
            molecule_file = user_input.get_file(prompt, standards.PDBDirectory, fetch_pdb=True)
            try:
                # Try to load the molecule.
                molecule = molecules.Molecule(name, molecule_file, self, exclude_hetero=self.args.exclude_hetero)
                molecule.get_structure()
            except Exception as error:
                # If that fails, say so and ask for another file name.
                disp("There was an error with the PDB import:\n{}".format(error.message))
            else:
                do = False
        return molecule

    def select_epitope_residues(self, chains):
        """ Ask the user to specify which residues are part of the antigen's target epitope. """
        self.epitope_residue_ids = dict()
        for chain in chains:
            # Ask for residues in each chain separately.
            self.epitope_residue_ids[chain.get_id()] = user_input.select_from_list("Please select the epitope residues from chain {}: ".format(chain.get_id()), [residue.get_id() for residue in chain], 1, None, map(molecules.residue_code, chain))

    def write_antigen_input_residues(self, entire_input_structure):
        """ Write a file containing only the atoms in the antigen.
        This function is called at the beginning of the Experiment to remove atoms in the input file that aren't in the antigen.
        """
        self.antigen_input_name = "antigen_input_file"
        # Make a file for the antigen atoms called antigen_input.pdb.
        self.antigen_input_file = os.path.join(self.structure_directory, os.path.basename("antigen_input.pdb"))
        # Make a molecule for the antigen. 
        self.antigen_input_molecule = molecules.Molecule(self.antigen_input_name, self.antigen_input_file, self)
        # Make a BioPython selector that will select only atoms in the antigen.
        selector = molecules.SelectChains(self.antigen_chain_ids)
        # Give the molecule all of the atoms and have it write only those in the antigen (using the selector).
        self.antigen_input_molecule.write_pdb(entire_input_structure, selector=selector)

    def safe_rmtree(self, directory):
        """ Safely remove a directory tree. """
        if not self.args.keeptemp:
            # For safety (i.e. no 'rm -rf /' accidents), ensure that the directory to be removed is in the Experiment's directory.
            if standards.is_subdirectory(directory, self.directory):
                try:
                    standards.safe_rmtree(directory)
                except:
                    #FIXME: This could potentially cause silent failures when the directory should be removed but is not.
                    pass
            else:
                raise OSError("{} cannot remove directory trees outside of {}".format(self.name, self.directory))

    def purge_temp(self):
        """ Remove the temporary directory. """
        standards.safe_rmtree(self.get_temp())

    def completed(self, args):
        """ Print a message saying that the Experiment has finished. """
        disp("{} has finished running. Please view the results in {}".format(self.name, self.directory))


class InteractionExperiment(Experiment):
    """ This Experiment simply calculates the interaction energy between two groups of chains, performs a relaxation, then calculates the energy again. """
    def __init__(self, args=None):
        Experiment.__init__(self, args)
        # FIXME: Allow the user to specify settings other than the defaults.
        self.topology_files = [standards.DefaultTopologyFile]
        self.parameter_files = [standards.DefaultParameterFile]
        self.solvation_files = [standards.DefaultSolvationFile]
        self.charmm_energy_terms = standards.DefaultCharmmEnergyTerms
        self.charmm_iterations = standards.DefaultCharmmIterations
        self.walltime = standards.DefaultWalltime
        # Specify the file from which to obtain the interacting molecules.
        self.molecule = self.get_molecule("input", "Please enter the name of the molecule for {}: ".format(self.name))
        self.molecule_file = self.molecule.file
        # Select chains for groups 1 and 2. The interaction energy between these groups will be calculated..
        self.chain_ids = [chain.get_id() for chain in self.molecule.get_chains()]
        chains_left = list(self.chain_ids)
        self.chain_groups = list()
        for i in [1, 2]:
            group = user_input.select_from_list("Please select the chain(s) for group {}: ".format(i), chains_left, min_number=1, max_number=(len(chains_left) + i - 2))
            for chain_id in group:
                chains_left.pop(chains_left.index(chain_id))
            self.chain_groups.append(group)
        # Initialize the status to 0 and submit the first script to the queue.
        self.status = 0
        self.save_and_submit()

    def get_task(self, status):
        """ Return the task (a function) and the purpose of the task (a string). There is just one task for an InteractionExperiment. """
        return [(self.relaxation, "Relaxation")][status]

    def relaxation(self, args):
        """ Relax the structure of the molecules before performing the interaction energy calculation. """
        self.energy_fields = ["Condition", "Complex", "Group 1", "Group 2", "Interaction"]
        garbage = list()
        # Calculate interaction energy before relaxation.
        group1, group2 = self.molecule.disassemble_into_chains(self.chain_groups)
        with charmm.InteractionEnergy(self, [group1, group2], solvation=False) as e:
            # Store the energy of the molecules separately and in complex.
            self.energies_unrelaxed = {"Condition": "Unrelaxed"}
            self.energies_unrelaxed["Group 1"] = e.energy1
            self.energies_unrelaxed["Group 2"] = e.energy2
            self.energies_unrelaxed["Complex"] = e.energy12
            self.energies_unrelaxed["Interaction"] = e.energy
            if self.benchmarking:
                self.add_drive_usage()
        # Relax the complex.
        self.relaxed_file = os.path.join(self.structure_directory, "relaxed.pdb")
        self.molecule.relax(relaxed_file=self.relaxed_file, in_place=True)
        # Calculate interaction energy after relaxation.
        group1, group2 = self.molecule.disassemble_into_chains(self.chain_groups)
        garbage.extend([group1.file, group2.file])
        with charmm.InteractionEnergy(self, [group1, group2], solvation=True) as e:
            # Store the energy of the molecules separately and in complex.
            self.energies_relaxed = {"Condition": "Relaxed"}
            self.energies_relaxed["Group 1"] = e.energy1
            self.energies_relaxed["Group 2"] = e.energy2
            self.energies_relaxed["Complex"] = e.energy12
            self.energies_relaxed["Interaction"] = e.energy
            if self.benchmarking:
                self.add_drive_usage()
        # Remove all temporary files.
        for fn in garbage:
            try:
                os.remove(fn)
            except OSError:
                pass
        # Write the report.
        self.create_report()
        # Remove the temporary file tree.
        self.safe_rmtree(self.get_temp())
    
    def create_report(self):
        """ Create a report for the InteractionExperiment. """
        summary = [
            ("Experiment name", self.name),
            ("Molecule input file", self.molecule_file),
            ("Group 1 chains", ", ".join(self.chain_groups[0])),
            ("Group 2 chains", "; ".join(self.chain_groups[1]))
        ]
        self.summary_file = os.path.join(self.directory, "Summary.txt")
        with open(self.summary_file, "w") as f:
            f.write("\n".join(["{}: {}".format(field, value) for field, value in summary]))
        self.results_file = os.path.join(self.directory, "Results.csv")
        with open(self.results_file, "w") as f:
            writer = csv.DictWriter(f, self.energy_fields)
            writer.writeheader()
            writer.writerow(self.energies_unrelaxed)
            writer.writerow(self.energies_relaxed)


class AntibodyInteractionExperiment(InteractionExperiment):
    """ This Experiment constructs an antibody-antigen complex given an antigen, an antigen position, and six MAPs parts; then it calculates the interaction energy, relaxes the complex, and recalculates the energy. """
    def __init__(self, args=None):
        Experiment.__init__(self, args)
        # FIXME: allow the user to specify settings other than the defaults.
        self.topology_files = [standards.DefaultTopologyFile]
        self.parameter_files = [standards.DefaultParameterFile]
        self.solvation_files = [standards.DefaultSolvationFile]
        self.charmm_energy_terms = standards.DefaultCharmmEnergyTerms
        self.charmm_iterations = standards.DefaultCharmmIterations
        self.walltime = standards.DefaultWalltime
        self.antigen = self.get_molecule("antigen_input", "Please enter the file of the antigen for {}: ".format(self.name))
        self.antigen_input_file = self.antigen.file
        # Select antigen chains.
        antigen_chains_all = {chain.get_id(): chain for chain in self.antigen.get_chains()}
        antigen_chains_selected = user_input.select_from_list("Please select the chain(s) for the antigen: ", antigen_chains_all, min_number=1)
        self.antigen_chain_ids = [chain.get_id() for chain in antigen_chains_selected]
        self.select_epitope_residues(antigen_chains_selected)
        # Create the antigen-antibody complex.
        self.antigen_position = [self.ask(pos, number=True) for pos in standards.PositionOrder]
        do = True
        while do:
            try:
                # Ask for the MAPs parts.
                parts = [self.ask("{} CDR".format(cdr)) for cdr in standards.MapsNamesakeCdrs]
                self.maps_parts = dict()
                for part in parts:
                    cdr, number = maps.split(part)
                    self.maps_parts[cdr] = number
                # Create a protoantibody from the MAPs parts.
                proto = self.get_proto_antibody()
            except Exception as e:
                disp(e.message)
            else:
                do = False
        # Initialize the status to 0 and submit the initial script to the queue.
        self.status = 0
        self.save_and_submit()

    def create_report(self):
        """ Create a report of the experiment results. """
        summary = [
            ("Experiment name", self.name),
            ("Antigen input file", self.antigen_input_file),
            ("Antigen chains", ", ".join(self.chain_groups[0])),
            ("Antibody chains", "; ".join(self.chain_groups[1])),
            ("Antigen position", ", ".join(map(str, self.antigen_position))),
            ("MAPs parts", ", ".join([maps.join(cdr, number) for cdr, number in self.maps_parts.items()]))
        ]
        self.summary_file = os.path.join(self.directory, "Summary.txt")
        with open(self.summary_file, "w") as f:
            f.write("\n".join(["{}: {}".format(field, value) for field, value in summary]))
        self.results_file = os.path.join(self.directory, "Results.csv")
        with open(self.results_file, "w") as f:
            writer = csv.DictWriter(f, self.energy_fields)
            writer.writeheader()
            writer.writerow(self.energies_unrelaxed)
            writer.writerow(self.energies_relaxed)

    def get_proto_antibody(self):
        """ Construct a ProtoAntibody object from the specified MAPs parts and the antigen position. """
        return molecules.ProtoAntibody(self.maps_parts, self.antigen_position, energy=None)

    def get_task(self, status):
        """ Return the task (a function) and the purpose (a str) of the experiment. """
        return [
            (self.antigen_relaxation, "Antigen relaxation and positioning"),
            (self.relaxation, "Relaxation")
        ][status]

    def antigen_relaxation(self, args):
        """ Relax the antigen structure before assembling the complex. """
        # Make a file of just the antigen structure.
        self.write_antigen_input_residues(self.antigen.get_structure())
        # Create and relax the antigen molecule.
        antigen_molecule = molecules.Molecule(self.antigen_input_name, self.antigen_input_file, self)
        self.antigen_relaxed_name = "relaxed_antigen"
        self.antigen_relaxed_file = os.path.join(self.structure_directory, "antigen_relaxed.pdb")
        antigen_molecule.relax(self.antigen_relaxed_file)
        # Rotate the epitope such that the z coordinate of its geometric center is minimized.
        self.epitope_zmin_file = self.antigen_relaxed_file
        antigen_molecule.position_antigen(0, 0, 0, 0, in_place=True)
        # Assemble and relax the antigen-antibody complex.
        self.create_complex()

    def create_complex(self):
        """ Assemble and relax the antigen-antibody complex. """
        ab_name = "antibody-antigen_complex"
        # Build a molecule from the proto-antibody.
        self.molecule = self.get_proto_antibody().to_molecule(ab_name, self.epitope_zmin_file, self)
        self.antibody_chain_ids = [chain.get_id() for chain in self.molecule.get_chains() if chain.get_id() not in self.antigen_chain_ids]
        # Partition the chains into antigen and antibody.
        self.chain_groups = list()
        self.chain_groups.append(self.antigen_chain_ids)
        self.chain_groups.append(self.antibody_chain_ids)
        # Relax the complex and compute the interaction energy.
        self.relaxation(None)
        

class ResidueContactExperiment(Experiment):
    """ Identify residues in the interface between two molecules. """
    def __init__(self, args=None):
        Experiment.__init__(self, args)
        # The user inputs the name of the file containing the input structure.
        self.molecule = self.get_molecule("input", "Please enter the name of the molecule for {}: ".format(self.name))
        self.chain_ids = [chain.get_id() for chain in self.molecule.get_chains()]
        chains_left = list(self.chain_ids)
        self.chain_groups = list()
        self.group_numbers = [1, 2]
        # The user specifies which chains are in the two molecules in contact.
        for i in self.group_numbers:
            group = user_input.select_from_list("Please select the chain(s) for group {}: ".format(i), chains_left, min_number=1, max_number=(len(chains_left) + i - 2))
            # Remove each chain that the user selected from the list of unselected chains so that the user cannot select the same chain twice.
            for chain_id in group:
                chains_left.pop(chains_left.index(chain_id))
            self.chain_groups.append(group)
        # Create a dictionary of {chain: group_number for each chain}.
        g1, g2 = self.chain_groups
        group_from_chain = dict()
        for group_number, group in zip(self.group_numbers, self.chain_groups):
            for chain in group:
                group_from_chain[chain] = group_number
        # The user specifies the maximum distance at which residues can be considered interacting.
        self.cutoff = self.ask("cutoff distance", number=True)
        # Find the inter-residue contacts.
        contacts = self.molecule.interchain_residue_contacts(g1, g2, self.cutoff)
        # Restructure the contact map so that it looks like {group_number: {chain: [residue, residue, ...], chain: ...}, group_number: ...}
        self.contacts_from_group = defaultdict(lambda: defaultdict(list))
        for chain_pair, chain_pair_contacts in contacts.iteritems():
            for chain in chain_pair:
                group = group_from_chain[chain]
                residues = self.contacts_from_group[group][chain]
                for contact in chain_pair_contacts:
                    residue = contact[chain]
                    if residue not in residues:
                        residues.append(residue)
        # Sort the residues by their ids, and then convert them into residue codes.
        for group, chains in self.contacts_from_group.iteritems():
            for chain, residues in chains.iteritems():
                residues.sort(key=Residue.Residue.get_id)
                for i in range(len(residues)):
                    residues[i] = molecules.residue_code(residues[i])
        # Write the report.
        self.create_report()

    def create_report(self):
        """ Create a report of the interacting residues. """
        summary = [
            ("Experiment name", self.name),
            ("Molecule input file", self.molecule.file),
            ("Group 1 chains", ", ".join(self.chain_groups[0])),
            ("Group 2 chains", "; ".join(self.chain_groups[1])),
            ("Cutoff distnce", self.cutoff)
        ]
        self.summary_file = os.path.join(self.directory, "Summary.txt")
        with open(self.summary_file, "w") as f:
            f.write("\n".join(["{}: {}".format(field, value) for field, value in summary]))
        self.results_file = os.path.join(self.directory, "Results.csv")
        with open(self.results_file, "w") as f:
            fields = ["Group", "Chain", "Residues"]
            writer = csv.DictWriter(f, fields)
            writer.writeheader()
            for group in self.group_numbers:
                for chain in self.chain_ids:
                    if chain in self.contacts_from_group[group]:
                        info = {
                            "Group": str(group),
                            "Chain": chain,
                            "Residues": ",".join(self.contacts_from_group[group][chain])
                        }
                        writer.writerow(info)


class TransformMoleculeExperiment(Experiment):
    """ This is a simple Experiment that just translates and rotates a molecule. """
    def __init__(self, args=None):
        Experiment.__init__(self, args)
        self.molecule = self.get_molecule("input", "Please enter the name of the molecule for {}: ".format(self.name))
        # Get the translation vector.
        do = True
        while do:
            try:
                translation_vector = np.array(map(float, self.ask("translation vector").split()))
            except ValueError:
                pass
            else:
                do = False
        # Get the rotation angle.
        rotation_degrees = self.ask("rotation angle (in degrees)", number=True)
        rotation_degrees -= np.floor(rotation_degrees / 360.0)
        # Get the rotation axis.
        do = not np.isclose(rotation_degrees, 0.0)
        rotation_axis = None
        while do:
            try:
                rotation_axis = np.array(map(float, self.ask("rotation axis").split()))
            except ValueError:
                pass
            else:
                do = False
        # Perform the translation.
        self.output_file = os.path.join(self.structure_directory, self.ask("output molecule", valid_path=True))
        self.molecule.translate(translation_vector, file_name=self.output_file, in_place=True)
        # Perform the rotation.
        if rotation_axis is not None:
            self.molecule.rotate(standards.rotate_axis_angle(rotation_axis, rotation_degrees, degrees=True), in_place=True)
        # The rotated molecules are in the experiment's structures directory. There is no Results directory.
        self.completed(None)


class OptmavenExperiment(Experiment):
    """ OptmavenExperiment designs de novo humanized antibody variable domains targeting a user-specified epitope. """
    def __init__(self, args=None):
        Experiment.__init__(self, args)
        # Ask whether to customize advanced settings.
        if user_input.get_yn("Would you like to customize Optmaven settings? "):
            self.grid_x = user_input.get_range("Antigen positioning grid, X axis")
            self.grid_y = user_input.get_range("Antigen positioning grid, Y axis")
            self.grid_z = user_input.get_range("Antigen positioning grid, Z axis")
            self.grid_zAngle = user_input.get_range("Antigen positioning grid, Z angle", 0, 360)
            self.clash_cutoff = user_input.get_number("Antigen positioning clash cutoff (Angstroms): ", 0, None, float)
            self.topology_files = user_input.get_files("Topology file(s): ", standards.InputsDirectory, 1, None)
            self.parameter_files = user_input.get_files("Parameter file(s): ", standards.InputsDirectory, 1, None)
            self.solvation_files = user_input.get_files("Solvation file(s): ", standards.InputsDirectory, None, None)
            self.charmm_energy_terms = user_input.select_from_list("CHARMM energy term(s): ", standards.DefaultCharmmEnergyTerms, 1, None)
            self.charmm_iterations = user_input.get_number("CHARMM iterations: ", 1, 20000, int)
            self.walltime = user_input.get_number("Walltime (seconds): ", 1, 35999, int)
            self.batch_size = user_input.get_number("Batch size: ", 1, None, int)
            self.number_of_designs = user_input.get_number("Number of designs output: ", 1, None, int)
        else:
            # If not, use the default settings.
            self.grid_x = standards.DefaultOptmavenGrid_x
            self.grid_y = standards.DefaultOptmavenGrid_y
            self.grid_z = standards.DefaultOptmavenGrid_z 
            self.grid_zAngle = standards.DefaultOptmavenGrid_zAngle
            self.clash_cutoff = standards.DefaultClashCutoff
            self.topology_files = [standards.DefaultTopologyFile]
            self.parameter_files = [standards.DefaultParameterFile]
            self.solvation_files = [standards.DefaultSolvationFile]
            self.charmm_energy_terms = standards.DefaultCharmmEnergyTerms
            self.charmm_iterations = standards.DefaultCharmmIterations
            self.walltime = standards.DefaultWalltime
            self.batch_size = standards.DefaultBatchSize
            self.number_of_designs = standards.DefaultNumberOfDesigns
        # Define the molecule containing the antigen.
        entire_input_model = self.get_entire_input_model()
        # Select the antigen chains.
        antigen_input_chains = self.select_antigen_input_chains(entire_input_model)
        self.antigen_chain_ids = [chain.get_id() for chain in antigen_input_chains]
        # Specify the epitope.
        self.select_epitope_residues(antigen_input_chains)
        # Write a file containing only the antigen molecule.
        self.write_antigen_input_residues(entire_input_model)
        # Initialize the status to 0, report the directory where the results will be, and submit the initial script to the queue.
        self.status = 0
        self.report_directory()
        self.save_and_submit()

    def change_status(self, new_status=None):
        """ Change the expriment to a given status. If no status is given, then increment the status. """
        if new_status is None:
            new_status = self.status + 1
        self.status = new_status
        # Save the experiment so that the new status is visible to all scripts.
        self.save()

    def get_task(self, status):
        """ Return the task (a function) and the status (a str) of the experiment, based on the current status. """
        return [
            (self.relax_antigen, "Antigen relaxation, positioning, and grid search"),
            (self.maps_energy_batch, "MAPs energy calculations"),
            (self.collect_maps_energies, "Prepare for MILP design"),
            (self.select_parts_batch, "MILP design"),
            (self.select_antibodies, "Select designs"),
            (self.relax_complexes_batch, "Relaxing designs"),
            (self.create_report, "Creating report"),
            (self.completed, "Completed")
        ][status]

    def create_report(self, args):
        """ Create a report of the results, including a summary and results file. """
        # Define the information to put in the Summary.txt file.
        summary = [
            ("Experiment name", self.name),
            ("Topology files", ", ".join(self.topology_files)),
            ("Parameter files", ", ".join(self.parameter_files)),
            ("Solvation files", ", ".join(self.solvation_files)),
            ("CHARMM energy terms", ", ".join(self.charmm_energy_terms)),
            ("CHARMM iterations", str(self.charmm_iterations)),
            ("Antigen input file", self.entire_input_file),
            ("Antigen chains", ", ".join(self.antigen_chain_ids)),
            ("Epitope residues", "; ".join(["{}: {}".format(chain, ", ".join(["".join(map(str, resid)) for resid in residues])) for chain, residues in self.epitope_residue_ids.items()])),
            (standards.zAngleLabel, ", ".join(map(str, self.grid_zAngle))),
            (standards.xLabel, ", ".join(map(str, self.grid_x))),
            (standards.yLabel, ", ".join(map(str, self.grid_y))),
            (standards.zLabel, ", ".join(map(str, self.grid_z))),
            ("Clash cutoff", self.clash_cutoff),
            ("Total positions", len(self.positions)),
            ("Selected designs", self.select_number)
        ]
        # Write the summary.
        self.summary_report = os.path.join(self.directory, "Summary.txt")
        with open(self.summary_report, "w") as f:
            f.write("\n".join(["{}: {}".format(field, value) for field, value in summary]))
        # Define a list to store results information.
        result_report_info = list()
        # Define a list of fields of results to report.
        result_report_fields = list()
        # Loop through all of the results (i.e. designs of antibodies).
        for index in range(self.select_number):
            # The results are stored as pickled OptmavenResult objects. Load them.
            with open(self.results_pickle_file(index)) as f:
                result = pkl.load(f)
            # Copy the information from the OptmavenResult into the result_info list.
            result_info = result.output()
            result_report_info.append(result_info)
        # Add any as-yet undefined fields in the results to the list of fields.
        for result in result_report_info:
            for field in result:
                if field not in result_report_fields:
                    result_report_fields.append(field)
        # Write the results.
        self.result_report = os.path.join(self.directory, "Results.csv")
        with open(self.result_report, "w") as f:
            writer = csv.DictWriter(f, fieldnames=result_report_fields)
            writer.writeheader()
            for result in result_report_info:
                writer.writerow(result)
        # Output benchmarking results, if applicable.
        if self.benchmarking:
            self.benchmarking_file = os.path.join(self.directory, "Benchmarking.csv")
            with open(self.benchmarking_file, "w") as f:
                # Define a writer object to do the writing.
                writer = csv.DictWriter(f, fieldnames=standards.BenchmarkingFields)
                # Initialize the total for each field to 0.
                totals = {field: 0.0 for field in standards.UnixTimeCodes.keys() + ["Drive Usage"]}
                totals["Type"] = "Total"
                writer.writeheader()
                # List the files storing benchmarking results.
                benchmark_files = [os.path.join(self.benchmark_directory, fn) for fn in os.listdir(self.benchmark_directory)]
                benchmarks = list()
                # Load the benchmarking results.
                for fn in benchmark_files:
                    # Each benchmarking result is stored as a pickled benchmarking Task object.
                    with open(fn) as f:
                        benchmarks.append(pkl.load(f))
                # Sort by time stamp so that the results appear in the order in which the tasks were performed.
                benchmarks.sort(key=benchmarking.Task.get_time_stamp)
                tried_current_time_file = False
                # Loop through the benchmarking Task objects to load their stored information.
                for benchmark in benchmarks:
                    try:
                        # Try to load the information.
                        d = benchmark.to_dict()
                    except benchmarking.BlankTimeFileError:
                        # Expect to find exactly one blank time file: the one measuring the current process.
                        if tried_current_time_file:
                            raise ValueError("Found two blank time files.")
                        tried_current_time_file = True
                    else:
                        # If the information was loaded successfully, then write it to the results file.
                        writer.writerow(d)
                        if isinstance(benchmark, benchmarking.Time):
                            # Add the current times to the total time for each time code (real, sys, and user), to keep track of the times.
                            for field in standards.UnixTimeCodes:
                                totals[field] += d[field]
                        elif isinstance(benchmark, benchmarking.DriveUsage):
                            # The Drive Usage total is the maximum drive usage, not the sum of all drive usages.
                            totals["Drive Usage"] = max(d["Drive Usage"], totals["Drive Usage"])
                # Remove the temporary directory.
                self.safe_rmtree(self.temp)
                # Record the drive usage after all temporary files have been removed.
                du = self.add_drive_usage(_file=False).to_dict()
                writer.writerow(du)
                totals["Drive Usage"] = max(du["Drive Usage"], totals["Drive Usage"])
                writer.writerow(totals)
        else:
            # If no benchmarking results must be output, then remove the temporary directory.
            self.safe_rmtree(self.temp) 

    def relax_antigen(self, args):
        """ Relax the antigen structure and then rotate it to position the epitope. """
        # Make a molecule for the antigen.
        antigen_molecule = molecules.Molecule(self.antigen_input_name, self.antigen_input_file, self)
        self.antigen_relaxed_name = "relaxed_antigen"
        self.antigen_relaxed_file = os.path.join(self.structure_directory, "antigen_relaxed.pdb")
        # Relax the antigen.
        antigen_molecule.relax(self.antigen_relaxed_file)
        # Position the epitope.
        self.minimize_epitope_z_coordinates(antigen_molecule)

    def minimize_epitope_z_coordinates(self, antigen_molecule, grid_search=True):
        """ Rotate the antigen around its geometric center so as to minimize the z coordinate of the epitope's geometric center. """
        self.epitope_zmin_file = self.antigen_relaxed_file
        # Rotate the epitope by calling the position_antigen function, which performs the rotation automatically.
        # position_antigen also moves the epitope's geometric center to the origin..
        antigen_molecule.position_antigen(0, 0, 0, 0, in_place=True)
        # Position the antigen using a grid search.
        if grid_search:
            self.grid_search()
        
    def grid_search(self):
        """ Search over the user-defined antigen-binding site and output all positions in which the antigen does not clash sterically with the antibody framework. """
        # Define a file to record the positions.
        self.positions_file = os.path.join(self.temp, "positions.dat")
        # Perform the grid search.
        with klaus.GridSearch(self) as x:
            if self.benchmarking:
                # If benchmarking, record the drive usage immediately after the grid search.
                self.add_drive_usage()
        # Increment the status.
        self.change_status()
        # Run the MAPs energy calculations.
        self.maps_energies_all(None)

    def get_maps_part_energy_directory_finished(self, part):
        return os.path.join(self.maps_energies_directory, part)

    def get_maps_part_energy_file_finished(self, part):
        return os.path.join(self.get_maps_part_energy_directory_finished(part), "{}_energies.dat".format(part))

    def maps_energies_all(self, args):
        """ Calculate the interacton energy between the antigen and all MAPs parts. """
        # Make a directory in which to store the interaction energy results.
        self.maps_energies_directory = os.path.join(self.temp, "maps_energies")
        try:
            os.mkdir(self.maps_energies_directory)
        except OSError:
            pass
        # Separate scripts will be created and submitted to actually calculate the interaction energies.
        # The experiment must first be saved so that these scripts have access to the following from the experiment's pickle file:
        # 1. the new status (set in self.grid_search) that indicates that MAPs interaction energies are to be calculated
        # 2. the path to the self.maps_energies_directory in which to store the energies.
        self.save()
        # The calculations are submitted as jobs, where each job calculates the interaction energy between one MAPs part and the antigen in all positions.
        # Thus there is one job per part.
        # jobs is a dict of {MAPs part: the file in which the interaction energies for that part will be saved.
        jobs = {part: self.get_maps_part_energy_file_finished(part) for part in maps.parts}
        # Submit all of the jobs using a PbsBatchSubmitter.
        self.submit(jobs=jobs)

    def maps_energy_batch(self, args):
        """ Calculate the interacton energy between the antigen and a batch of MAPs parts.
        The PbsBatchSubmitter (called from maps_energies_all) creates a number of so-called batches, which are sets of jobs.
        Each batch is run using its own script, which has the command 'path/to/python path/to/experiment.py path/to/experiment_name.pickle path/to/parts_file'.
        When python runs experiment.py (see if __name__ == '__main__' at the bottom), experiment.py unpickles the experiment (argument 1).
        The experiment.run function is then called; since the status is currently 'MAPs energy calculations,' maps_energy_batch is called.
        The arguments are passed to maps_energy_batch, which opens the parts_file (argument 2).
        The parts file is a text file listing the MAPs parts that should be run in this batch.
        """
        # Open the parts file and get a list of the MAPs parts to use.
        parts_file = args[2]
        with open(parts_file) as f:
            parts = f.read().split()
        # Calculate the energies for each part.
        for part in parts:
            self.maps_energy_single(part)
   
    def maps_energy_single(self, part):
        """ Calculate the interaction energy between one MAPs part and the antigen in all positions. """
        with klaus.MapsEnergies(self, part, self.get_maps_part_energy_file_finished(part)) as energies:
            if self.benchmarking:
                # Record the drive usage after the calculation.
                self.add_drive_usage()
    
    def collect_maps_energies(self, args):
        """ Load the MAPs interaction energies into a dict, then initiate the select parts step. """
        # Define a dict to store the maps energies.
        maps_energies = defaultdict(list)
        # Loop through all of the MAPs parts.
        for part in maps.parts:
            # Generate the name of the file containing the energies for the part.
            _file = self.get_maps_part_energy_file_finished(part)
            with open(_file) as f:
                # Read the antigen position and energy on each line.
                for line in f:
                    try:
                        # Each line should have five numbers: zAngle, x, y, z, and energy, in that order. 
                        zAngle, x, y, z, energy = map(float, line.split())
                    except ValueError:
                        raise ValueError("Cannot parse line in MAPs energy file {}:\n{}".format(_file, line))
                    # Convert the antigen position into a tuple that can be used as a hash.
                    position = (zAngle, x, y, z)
                    # Save the energies for each position as a list of 3-item lists.
                    # Each 3-item list contains the MAPs part CDR (e.g. LV), number (e.g. 3), and interaction energy at the corresponding position.
                    # Example: maps_energies = {(0.0, 0.0, 0.0, 0.0): [['HV', 1, -3.4], ['HV', 2, 10.3], ...], (60.0, 2.5, -5.0, 16.25): [['KCDR3', 115, -4.8], ...], ...}
                    # Grouping by position (rather than part) is important because each position will be used to generate one design.
                    maps_energies[position].append([maps.get_part_cdr(part), maps.get_part_number(part), energy])
        # Change the status to indicate that the experiment is in the 'Prepare for MILP design' step.
        self.change_status()
        # Prepare the scripts for selecting all parts (i.e. generating all designs).
        self.select_parts_all(maps_energies)

    def get_select_parts_directory(self, index):
        return os.path.join(self.select_parts_directory, "position_{}".format(index))

    def get_select_parts_energy_file(self, index):
        return os.path.join(self.get_select_parts_directory(index), "energies.pickle")

    def get_select_parts_file_finished(self, index):
        return os.path.join(self.get_select_parts_directory(index), "parts.pickle")

    def select_parts_all(self, maps_energies):
        """ Prepare the scripts for the part selection (i.e. design) step. """
        # Make a directory to store the designs.
        self.select_parts_directory = os.path.join(self.temp, "select_parts")
        try:
            os.mkdir(self.select_parts_directory)
        except OSError:
            pass
        # Index all of the positions so that they can each be described with an integer.
        # For example, self.positions = {0: (120.0, 5.0, -2.5, 5.0), 1: (300.0, 7.5, 5.0, 6.25), ...}
        # Indexing makes writing the job scripts easier because each position can be represented as one integer instead of a tuple of four floats.
        self.positions = dict(enumerate(maps_energies))
        # Pickle the MAPs energies so that the select parts scripts can use them.
        for index, position in self.positions.iteritems():
            # Pickle the energies for each position separately.
            position_energies = maps_energies[position]
            try:
                os.mkdir(self.get_select_parts_directory(index))
            except OSError:
                pass
            with open(self.get_select_parts_energy_file(index), "w") as f:
                pkl.dump(position_energies, f)
        if self.benchmarking:
            # Record the drive usage after saving all of these pickle files.
            self.add_drive_usage()
        # The experiment must be saved before submitting the jobs for the same reason it needed to be saved in maps_energies_all.
        self.save()
        # Remove maps energies directory: all energies are saved in pickle files.
        # NOTE: maps_energies was passed as an argument to select_parts_all rather than made an attribute of the experiment.
        # MAPs energies is a very large dict; if it were an attribute, then the experiment's pickle file would be enormous.
        # Thus, to avoid this storage burden, maps_energies is passed as an argument.
        # The energies are, however, saved as pickle files: one per position.
        # It would also be possible to save them all as a maps_energies attribute of the experiment.
        # However, each position needs to load just the maps energies for the position, not all of the energies.
        # Thus, each position would waste time unpickling the entire, enormous dict of maps energies, rather than the much smaller position-specific maps energies.
        self.safe_rmtree(self.maps_energies_directory)
        if self.benchmarking:
            # Record the drive usage after removing the maps energies directory.
            self.add_drive_usage()
        # Index the jobs in a manner analogous to the indexing in maps_energies_all.
        # It is much easier to index the positions than to use the actual position tuples because the index must be written into the filename, which is cumbersome with tuples.
        jobs = {index: self.get_select_parts_file_finished(index) for index in self.positions}
        # Submit the jobs.
        self.submit(jobs=jobs)

    def select_parts_batch(self, args):
        """ Select the parts for a batch of positions. This function works in a manner analogous to that of maps_energies_batch. """
        index_file = args[2]
        with open(index_file) as f:
            indexes = f.read().split()
        for index in indexes:
            self.select_parts_single(index)
        
    def select_parts_single(self, index, solution_cuts=None):
        """ Select the parts (generate the design) for one position. """
        index = int(index)
        position = self.positions[index]
        # Load the integer cuts
        clash_cuts = maps.get_integer_cuts()
        if solution_cuts is None:
            solution_cuts = list()
        with open(self.get_select_parts_energy_file(index)) as f:
            position_energies = pkl.load(f)
        # Load the MAPS energies
        # Using cplex to get the optimal solutions. 
        # The solution is a dictionary and the key and values are partname and number respectively
        # The energy is the total interaction energy between the selected parts and antigen
        solution, energy = maps.select_parts(position_energies, clash_cuts, solution_cuts)
        antibody = molecules.ProtoAntibody(solution, position, energy)
        with open(self.get_select_parts_file_finished(index), "w") as f:
            pkl.dump(antibody, f)
        if self.benchmarking:
            # Record the drive usage after generating this design.
            self.add_drive_usage()

    def select_antibodies(self, args):
        """ Cluster the antibodies based on their coordinates. """
        # Make a dict to store the designs.
        antibodies = defaultdict(list)
        # Loop through the design for each position.
        for design in map(self.get_select_parts_file_finished, self.positions):
            with open(design) as f:
                # Load the design, which is saved as a pickled ProtoAntibody object.
                antibody = pkl.load(f)
            # Bin the antibodies into the list for either kappa or lambda chain antibodies.
            # The kappa and lambda chain antibodies are clustered separately.
            antibodies[antibody.light_chain].append(antibody)
        clusters = dict()
        # For each chain (kappa and lambda):
        for chain, chain_abs in antibodies.iteritems():
            # Cluster the antibodies in the chain (listed in chain_abs).
            chain_clusters, mses = kmeans.optimal_kmeans(chain_abs)
            # Sort the designs in each cluster in order of increasing interaction energy.
            chain_clusters = [sorted(cluster) for cluster in chain_clusters]
            # Assign to each cluster an index, based on the minimum interaction energy of the cluster.
            clusters[chain] = dict(enumerate(sorted(chain_clusters)))
            # Write a file of the mse for each k during clustering.
            with open(os.path.join(self.get_temp(), "clustering_{}.csv".format(chain)), "w") as f:
                f.write("\n".join(["k,mse"] + ["{},{}".format(k, mse) for k, mse in mses.items()]))
        # Select the best antibodies from the clusters.
        self.highest_ranked_designs = list()
        # If there are more designs than the user wants, select the number that the user wants.
        # Otherwise, select all designs.
        self.select_number = min(self.number_of_designs, len(self.positions))
        while len(self.highest_ranked_designs) < self.select_number:
            # Collect in cluster_heads the best as-yet unselected design from each cluster that has at least one unselected design.
            cluster_heads = list()
            # Select separately from kappa and lambda clusters.
            for chain, chain_clusters in clusters.iteritems():
                # Loop through all the clusters for the light chain type, from the cluster with the lowest minimum interaction energy to that with the highest minimum interaction energy.
                for index in sorted(chain_clusters, key=lambda x: chain_clusters[x]):
                    cluster = chain_clusters[index]
                    if len(cluster) > 0:
                        # If a cluster has an unselected design, put the best design from the cluster in the 
                        antibody = cluster.pop(0)
                        antibody.set_cluster(index)
                        cluster_heads.append(antibody)
            # Add these designs to the list of best designs.
            while len(self.highest_ranked_designs) < self.select_number and len(cluster_heads) > 0:
                self.highest_ranked_designs.append(cluster_heads.pop(0))
        # Make a directory in which to store the unrelaxed antibody-antigen complexes for the best designs..
        self.unrelaxed_complex_directory = os.path.join(self.get_temp(), "unrelaxed_complexes")
        try:
            os.mkdir(self.unrelaxed_complex_directory)
        except OSError:
            pass
        # Define the names of the files of the unrelaxed complexes.
        self.unrelaxed_complex_files = [os.path.join(self.unrelaxed_complex_directory, "complex_{}.pdb".format(i)) for i in range(self.select_number)]
        # Create molecule objects for the unrelaxed complexes.
        [ab.to_molecule("unrelaxed", _file, self) for ab, _file in zip(self.highest_ranked_designs, self.unrelaxed_complex_files)]
        if self.benchmarking:
            # Record the drive usage after creating the molecules.
            self.add_drive_usage()
        # Remove the select parts directory; it is no longer needed.
        self.safe_rmtree(self.select_parts_directory)
        if self.benchmarking:
            # Record the drive usage after removing all of the pickled MAPs energies and unselected designs.
            self.add_drive_usage()
        # Change the status to 'Relaxing designs'
        self.change_status()
        # Relax all of the complexes.
        self.relax_complexes_all()

    def relaxed_complex_directory(self):
        return os.path.join(self.directory, "antigen-antibody-complexes")

    def relaxed_complex_file(self, index):
        return os.path.join(self.relaxed_complex_directory(), "complex_{}.pdb".format(index))

    def results_directory(self, index):
        return os.path.join(self.relaxed_complex_directory(), "Result_{}".format(index))

    def results_pickle_directory(self):
        return os.path.join(self.get_temp(), "results")

    def result_name(self, index):
        return "Result_{}".format(index)

    def results_pickle_file(self, index):
        return os.path.join(self.results_pickle_directory(), "result_{}.pickle".format(index))

    def relax_complexes_all(self):
        """ Relax all of the antigen-antibody complexes of the best designs. """
        # Make a directory into which the relaxed structure files will go.
        try:
            os.mkdir(self.relaxed_complex_directory())
        except OSError:
            pass
        # Make a directory for pickle files of the results.
        try:
            os.mkdir(self.results_pickle_directory())
        except OSError:
            pass
        # Create a job for each design to be relaxed.
        jobs = {index: self.results_pickle_file(index) for index in range(self.select_number)}
        # Save the experiment.
        self.save()
        # Submit the jobs.
        self.submit(jobs=jobs)

    def relax_complexes_batch(self, args):
        """ Relax a batch of antigen-antibody complexes. This function works much like maps_energies_batch. """
        # Get the indexes of the complexes to relax.
        index_file = args[2]
        with open(index_file) as f:
            indexes = map(int, f.read().split())
        # Relax each complex.
        for index in indexes:
            self.relax_complex_single(index)

    def relax_complex_single(self, index):
        """ Relax a single antigen-antibody complex. """
        garbage = list()
        _complex = molecules.AntibodyAntigenComplex("relaxed", self.unrelaxed_complex_files[index], self)
        # Calculate interaction energy before relaxation.
        antigen, antibody = _complex.disassemble_into_chains([_complex.antigen_chains, _complex.antibody_chains])
        garbage.extend([antigen.file, antibody.file])
        with charmm.InteractionEnergy(self, [antigen, antibody], solvation=False) as e:
            unrelaxed_energy = e.energy
            if self.benchmarking:
                self.add_drive_usage()
        # Relax the complex.
        _complex.relax(relaxed_file=self.relaxed_complex_file(index), in_place=True)
        # Calculate interaction energy after relaxation.
        antigen, antibody = _complex.disassemble_into_chains([_complex.antigen_chains, _complex.antibody_chains])
        garbage.extend([antigen.file, antibody.file])
        with charmm.InteractionEnergy(self, [antigen, antibody], solvation=True) as e:
            relaxed_energy = e.energy
            if self.benchmarking:
                self.add_drive_usage()
        # Create an OptmavenResult to store information about the relaxation (e.g. energy before and after relaxation).
        # This information will ultimately go into Results.csv.
        result = OptmavenResult(self.result_name(index), self.results_directory(index), _complex, self.highest_ranked_designs[index], unrelaxed_energy, relaxed_energy) 
        with open(self.results_pickle_file(index), "w") as f:
            pkl.dump(result, f)
        
    def get_entire_input_model(self):
        """ Ask the user to specify the file containing the antigen. """
        # Select the antigen file.
        self.entire_input_name = "entire_input_file"
        do = True
        while do:
            # Ask for the file name. If it doesn't exist but is a valid PDB ID, fetch the PDB.
            self.entire_input_file = user_input.get_file("Please name the file containing the antigen, or enter a PDB ID: ", standards.PDBDirectory, fetch_pdb=True)
            try:
                # Try to make a molecule from the file.
                entire_input_model = molecules.Molecule(self.entire_input_name, self.entire_input_file, self, exclude_hetero=self.args.exclude_hetero).get_model()
            except Exception as error:
                disp("There was an error with the PDB import:\n{}".format(error.message))
            else:
                do = False
        return entire_input_model

    def select_antigen_input_chains(self, entire_input_model):
        """ Specify which chains of the input molecule are part of the antigen. """
        chains = Selection.unfold_entities(entire_input_model, "C")
        chain_ids = [chain.get_id() for chain in chains]
        do = True
        while do:
            # Ask the user to select the chains.
            selected_chains = user_input.select_from_list("Please select the antigen chains: ", chains, 1, None, names=chain_ids)
            # When the antibody-antigen complexes are constructed, the antibody chains will be named H and L or K.
            # To prevent problems at this step, the antigen chain may not also be named one of these letters.
            # The user must manually rename the input molecule if the antigen contains a chain named H, K, or L.
            # FIXME: potentially include automatic renaming in the future.
            if any(chain.get_id() in standards.MapsChains for chain in selected_chains):
                disp("The following are not valid names for antigen chains: {}".format(", ".join(standards.MapsChains)))
            else:
                do = False
        return selected_chains


class OptmavenResult(object):
    """ Stores information about one result (a design) from an OptmavenExperiment. """
    def __init__(self, name, directory, molecule, proto_antibody, unrelaxed_energy, relaxed_energy):
        # Store the information.
        self.name = name  # name of the design
        self.directory = directory  # directory in which the design is located
        self.molecule = molecule  # molecule object for the design
        self.proto_antibody = proto_antibody  # proto-antibody object for the design
        self.unrelaxed_energy = unrelaxed_energy  # interaction energy of the complex before relaxation
        self.relaxed_energy = relaxed_energy  # interaction energy of the complex after relaxation

    def output(self):
        """ Output the information as a dict whose keys match the fields of Results.csv. """
        try:
            os.mkdir(self.directory)
        except OSError:
            pass
        info = OrderedDict()
        info["Result"] = self.name
        # Write the molecule.
        self.molecule_file = os.path.join(self.directory, "{}.pdb".format(self.name))
        shutil.copyfile(self.molecule.file, self.molecule_file)
        info["PDB file"] = self.molecule_file
        # Get the sequence and output it as a fasta.
        self.fasta_file = os.path.join(self.directory, "{}.fasta".format(self.name))
        records = list(SeqIO.parse(self.molecule_file, "pdb-atom"))
        SeqIO.write(records, self.fasta_file, "fasta")
        # Store the information about the result as a dict (info).
        for record in records:
            chain = record.id.split(":")[1]
            info[chain] = str(record.seq)
        for cdr, part in self.proto_antibody.get_namesake_parts().items():
            info[cdr] = part
        for dimension, coord in self.proto_antibody.get_labeled_position().items():
            info[dimension] = coord
        info["Cluster"] = self.proto_antibody.cluster
        info["MILP energy (kcal/mol)"] = self.proto_antibody.energy
        info["Unrelaxed energy (kcal/mol)"] = self.unrelaxed_energy
        info["Relaxed energy (kcal/mol)"] = self.relaxed_energy
        return info


class CreateAntigenAntibodyComplexExperiment(OptmavenExperiment):
    """ This Experiment creates an antigen-antibody complex without relaxing the structure. """
    def __init__(self, args):
        Experiment.__init__(self, args)
        # Ask the user to provide the input files and epitope residues.
        entire_input_model = OptmavenExperiment.get_entire_input_model(self)
        antigen_input_chains = OptmavenExperiment.select_antigen_input_chains(self, entire_input_model)
        self.antigen_chain_ids = [chain.get_id() for chain in antigen_input_chains]
        OptmavenExperiment.select_epitope_residues(self, antigen_input_chains)
        # Position the antigen.
        do = True
        while do:
            try:
                zAngle, x, y, z = map(float, self.ask("antigen position").split())
            except ValueError:
                pass
            else:
                do = False
        position = (zAngle, x, y, z)
        # Get the chain loci.
        heavy = user_input.select_one_from_list("Please specify the heavy chain: ", standards.MapsHeavyChains)
        light = user_input.select_one_from_list("Please specify the light chain: ", standards.MapsLightChains)
        # Get the MAPs parts.
        parts = dict()
        for cdr in standards.MapsCdrs:
            chain = maps.get_chain(cdr)
            if chain not in [heavy, light]:
                continue
            do = True
            while do:
                try:
                    number = int(self.ask(cdr))
                    part = maps.join(cdr, number)
                except ValueError as e:
                    disp(e.message)
                else:
                    do = False
            parts[cdr] = number
        energy = None
        output_file = user_input.get_file("Please specify the output file: ", self.directory, new_file=True)
        # Write the antigen input file.
        OptmavenExperiment.write_antigen_input_residues(self, entire_input_model)
        antigen_molecule = molecules.Molecule(self.antigen_input_name, self.antigen_input_file, self) 
        # Center and point the epitope downward.
        OptmavenExperiment.minimize_epitope_z_coordinates(self, antigen_molecule, grid_search=False)
        # Assemble the complex.
        self.proto_antibody = molecules.ProtoAntibody(parts, position, energy)
        self.proto_antibody.to_molecule("complex", output_file, self)
        clear()
        self.report_directory()


class ExperimentArgParser(argparse.ArgumentParser):
    """ Parse the command-line arguments and give these arguments to an Experiment. """
    def __init__(self):
        argparse.ArgumentParser.__init__(self)
        # Specify all of the available arguments and what they do.
        self.add_argument("--antigen", help="input antigen file")
        self.add_argument("--epitope", help="comma-separated list of epitope residues")
        self.add_argument("--exclude_hetero", help="specify whether to have [no] hetatm exclusion, [yes] exclude hetatms, or [ask] each time")
        benchtemp = self.add_mutually_exclusive_group()
        benchtemp.add_argument("--benchmark", help="turn benchmarking [on] or [off]")
        benchtemp.add_argument("--keeptemp", help="do dot remove maps energies or part files", action="store_true")

    def parse_args(self):
        """ Perform argument parsing. """
        args = argparse.ArgumentParser.parse_args(self)
        # If exclude_hetero is not specified, set it to standards.HetAtmAsk.
        if args.exclude_hetero is None:
            args.exclude_hetero = standards.HetAtmAsk
        # Ensure that exclude_hetero has a valid value.
        if args.exclude_hetero not in standards.HetAtmOptions:
            raise ValueError("Bad value for exclude_hetero: {}".format(args.exclude_hetero))
        return args


if __name__ == "__main__":
    """ If experiment.py is run directly, unpickle an run an Experiment.
    
    usage:
    path/to/python path/to/experiment.py path/to/experiment_name.pickle [arg1] [arg2] ...
    
    The experiment stored in experiment_name.pickle is unpickled, and then its run method is called.
    Any additional arguments are passed to run.
    """
    if len(sys.argv) < 2:
        raise OSError("Usage: path/to/python path/to/experiments.py path/to/experiment.pickle [arg1] [arg2] ...")
    # Unpickle the experiment.
    experiment_file = sys.argv[1]
    with open(experiment_file) as f:
        experiment = pkl.load(f)
    # Run the experiment with the arguments.
    experiment.run(sys.argv)
