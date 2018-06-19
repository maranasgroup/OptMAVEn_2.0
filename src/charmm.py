from collections import OrderedDict
import os
import shutil
import string
import subprocess
import sys
import tempfile

from Bio.PDB import Selection

import molecules
import standards


class CharmmProc(object):
    """ This class is the base class for managing all CHARMM procedures. """
    def __init__(self, experiment, molecule_list):
        self.experiment = experiment
        # One or multiple molecules can be involved in the CHARMM calculation.
        self.molecules = molecule_list
        # The residues of the molecules will be renumbered (starting at 1) before running CHARMM, or else CHARMM will crash.
        # We want to preserve the original residue numbering so that we can keep track of e.g. which residues are part of the specified epitope.
        # Thus, store all of the original residue numbers in the resnums OrderedDict.
        # resnums will be used to restore the original residue numbering after the CHARMM procedure.
        # NOTE: currently, this procedure is not equipped to handle residues with letters in their numbers, e.g. '101B'.
        self.resnums = OrderedDict()
        for molecule in self.molecules:
            molecule_resnums = OrderedDict()
            for model in molecule.get_structure():
                model_resnums = OrderedDict()
                for chain in model:
                    model_resnums[chain.get_id()] = [residue.get_id()[1] for residue in chain]
                molecule_resnums[model.get_id()] = model_resnums
            self.resnums.update(molecule_resnums)
        # Indicate that the CHARMM script has not yet run.
        self.has_run = False

    def begin_script(self):
        self.lines = introduction()

    def load_input_files(self):
        """ Load the Topology and Parameter input files in a CHARMM script. """
        # Loop through topology and parameter data.
        input_files = [
            {"title": "topology" , "files": self.experiment.topology_files,  "cmd": "rtf" },
            {"title": "parameter", "files": self.experiment.parameter_files, "cmd": "para"}
        ]
        for file_type in input_files:
            self.lines.append("! Load the {} file(s)".format(file_type["title"]))
            for count, file_ in enumerate(file_type["files"]):
                # Make a symbolic link in the current directory to the file.
                if not os.path.isfile(file_):
                    raise OSError("File does not exist: {}".format(file_))
                file_name = os.path.basename(file_).lower()
                symlink = os.path.join(self.directory, file_name)
                os.symlink(file_, symlink)
                self.lines.append("open read unit 10 form name {} card".format(file_name))
                # If there is more than one file, indicate that the subsequent files are being appended.
                append = " append" * int(count > 0)
                self.lines.append("read {} card unit 10{}".format(file_type["cmd"], append))
                self.lines.append("close unit 10")

    def load_molecules(self):
        """ Create text to load Molecules in a CHARMM script. """
        # Merge all of the molecules.
        merged_structure = molecules.merge(self.molecules)
        # Save the merged structure in charmm_format_file.
        charmm_format_name = "charmm_format"
        handle, charmm_format_file = tempfile.mkstemp(suffix=".pdb", dir=self.directory)
        merged_molecule = molecules.Molecule(charmm_format_name, charmm_format_file, self.experiment)
        # CHARMM needs all of the molecules to be in separate PDB files, so disassemble the merged molecule into separate files.
	chains, files = merged_molecule.disassemble_for_CHARMM(structure=merged_structure, directory=self.directory)
        # Remove the file of the merged structure, which is no longer needed.
        try:
            os.remove(charmm_format_file)
        except OSError:
            pass
        # Write the part of the script that loads the chains.
        self.chain_ids = list()
        for chain, _file in zip(chains, files):
            _id = chain.get_id()
            self.chain_ids.append(_id)
            segment = make_segment_name(_id)
            base_name = os.path.basename(_file)
            self.lines.extend(["! Load Chain {}".format(_id),
                "open read unit 10 form name {}".format(base_name),
                "read sequ pdb offi unit 10",
                "close unit 10",
                "gene {} setup".format(segment),
                "open read unit 10 form name {}".format(base_name),
                "read coor pdb unit 10",
                "close unit 10"])
        # Add any missing atoms to the chain.
        self.lines.extend(ic_fill())

    def calculate_energy(self, solvation=True):
        """ Create text to calculate the energy in a CHARMM script. """
        self.energy_file = os.path.join(self.directory, "energy.dat")
        self.lines.extend(["! Calculate the energy.",
            energy_line(self.experiment, solvation),
            "ener",
            "set tot ?ener",
            "open write card unit 10 name {}".format(os.path.basename(self.energy_file)),
            "write title unit 10",
            "*@tot",
            "close unit 10"])

    def relax(self, solvation):
        """ Create text to perform a structural relaxation in a CHARMM script. """
        self.lines.extend(["! Carry out an energy minimization",
            "nbon nbxm 5",
            energy_line(self.experiment, solvation),
            "mini abnr nstep {} nprint 50 -".format(self.experiment.charmm_iterations),
            "tolgrd 0.01 tolenr 0.0001 tolstp 0.00"])

    def output_molecules(self):
        """ Create text to output the molecules in a CHARMM script. """
        self.output_molecule_files = list()
        # NOTE: this function must be called AFTER self.load_molecules because self.load_molecules defines self.chain_ids.
        for _id in self.chain_ids:
            base_name = make_pdb_name(_id, "out")
            self.output_molecule_files.append(os.path.join(self.directory, base_name))
            segment = make_segment_name(_id)
            self.lines.extend(["open write unit 10 name {} card".format(base_name),
                "write coor sele segi {} end pdb unit 10 card".format(segment),
                "close unit 10"])

    def end_script(self):
        """ Write a STOP indicator at the end of a CHARMM script. """
        self.lines.append("STOP")

    def charmm(self):
        """ Run CHARMM using the script. """
        # Write the input script.
        self.script_file = os.path.join(self.directory, "charmm_input.txt")
        with open(self.script_file, "w") as f:
            f.write("\n".join(self.lines))
        # Create the output file.
        self.output_file = os.path.join(self.directory, "charmm_output.txt")
        # Simultaneously open the input and output scripts, and call CHARMM with subprocess.Popen.
        with open(self.script_file) as fi, open(self.output_file, "w") as fo:
            # This command securely directs the contents of the input script into CHARMM and directs the CHARMM output into the output file.
            proc = subprocess.Popen(standards.CharmmCommand, stdin=fi, stdout=fo)
            # Wait for CHARMM to finish.
            proc.wait()
        # Indicate that CHARMM has run.
        self.has_run = True

    def collect_garbage(self):
        """ Remove the directory associated with this CHARMM process. """
        self.experiment.safe_rmtree(self.directory)

    def __enter__(self):
        """ Automatically create a directory in which all of the actions of this process will occur. """
        self.directory = tempfile.mkdtemp(prefix="charmm_", dir=self.experiment.get_temp())
        self.previous_directory = os.getcwd()
        os.chdir(self.directory)

    def __exit__(self, exc_type, exc_value, traceback):
        """ Automatically delete the directory in which all of the actions of this process occured. """
        os.chdir(self.previous_directory)
        self.collect_garbage()


class Energy(CharmmProc):
    """ A CHARMM process that calculates the energy of molecules. """
    def __init__(self, experiment, molecule_list, solvation=True):
        CharmmProc.__init__(self, experiment, molecule_list)
        self.solvation = solvation

    def __enter__(self):
        # Specify the list of events that need to happen to calculate the energy.
        CharmmProc.__enter__(self)
        self.begin_script()
        self.load_input_files()
        self.load_molecules()
        self.calculate_energy(solvation=self.solvation)
        self.end_script()
        self.charmm()
        with open(self.energy_file) as f:
            try:
                self.energy = float(f.read())
            except ValueError:
                raise IOError("There was a problem with the CHARMM energy file {}".format(self.energy_file))
        return self


class InteractionEnergy(CharmmProc):
    """ A CHARMM process that calculates the interaction energy between two molecules. """
    def __init__(self, experiment, molecule_list, solvation=True):
        CharmmProc.__init__(self, experiment, molecule_list)
        # Make sure exactly two molecules have been specified.
        if len(self.molecules) == 2:
            self.mol1, self.mol2 = self.molecules
        else:
            raise ValueError("An interaction energy calculation needs exactly two molecules.")
        self.solvation = solvation

    def __enter__(self):
        # Interaction energy is calculated as the energy of the complex minus the energy of the components.
        # This isn't as efficient as it could be (i.e. using INTE) because all internal energies are computed, too.
        # However, this formulation allows calculation of solvation energy contribution to the interaction energy, which (to my knowledge) cannot be done using INTE.
        CharmmProc.__enter__(self)
        with Energy(self.experiment, [self.mol1], solvation=self.solvation) as e1:
            self.energy1 = e1.energy
        with Energy(self.experiment, [self.mol2], solvation=self.solvation) as e2:
            self.energy2 = e2.energy
        with Energy(self.experiment, self.molecules, solvation=self.solvation) as e12:
            self.energy12 = e12.energy
        self.energy = self.energy12 - (self.energy1 + self.energy2)
        return self


class Relaxation(CharmmProc):
    """ A CHARMM process that performs a structural relaxation. """
    def __init__(self, experiment, molecules, relaxed_file, ignore_solvation_initially=True):
        CharmmProc.__init__(self, experiment, molecules)
        self.relaxed_file = relaxed_file
        self.ignore_solvation_initially = ignore_solvation_initially

    def __enter__(self):
        # The events that need to happen during the relaxation.
        CharmmProc.__enter__(self)
        self.begin_script()
        self.load_input_files()
        self.load_molecules()
        # Relaxations with solvation can crash if there are steric clashes. Relax without solvation first to resolve any clashes, then relax again with solvation if needed.
        if self.ignore_solvation_initially or standards.CharmmSolvationTerm not in self.experiment.charmm_energy_terms:
            self.relax(solvation=False)
        if standards.CharmmSolvationTerm in self.experiment.charmm_energy_terms:
            self.relax(solvation=True)
        self.output_molecules()
        self.end_script()
        self.charmm()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        relaxed_name = "relaxed"
        relaxed_molecule = molecules.Molecule(relaxed_name, self.relaxed_file, self.experiment)
        relaxed_molecule.assemble_from_CHARMM(self.chain_ids, self.output_molecule_files, self.resnums)
        CharmmProc.__exit__(self, exc_type, exc_value, traceback)


def introduction(purpose=None):
    """Create a header to start a CHARMM script."""
    return ["wrnl -2",
        "!prnl -2",
        "bomb -2"]
        

def ic_fill():
    """Tell CHARMM to include missing Atoms."""
    return ["! Add missing Atoms and assign them coordinates",
        "ic fill preserve",
        "ic param",
        "ic build",
        "hbuild"]


def energy_line(experiment, solvation=True):
    """Generate a string saying how energy should be calculated in CHARMM."""
    terms = list(experiment.charmm_energy_terms)
    # If solvation is not being used, remove solvation from the terms list, if it exists.
    if not solvation:
        try:
            terms.pop(terms.index(standards.CharmmSolvationTerm))
        except ValueError:
            pass
    # If Generalized Born solvation is used, specify solvation parameters.
    if standards.CharmmSolvationTerm in terms:
        line = """GBMV BETA -20 EPSILON 80 DN 1.0 watr 1.4 GEOM -
    TOL 1e-8 BUFR 0.5 Mem 10 CUTA 20 HSX1 -0.125 HSX2 0.25 -
    ALFRQ 1 EMP 0.25 P4 0.0 P6 8.0 P3 0.70 ONX 1.9 OFFX 2.1 -
    WTYP 2 NPHI 38 SHIFT -0.102 SLOPE 0.9085 CORR 1
"""
    else:
        line = ""
    line += "skip all excl {}".format(" ".join(terms))
    return line


def make_segment_name(_id):
    """ Make a segment name from a chain ID. Simply make the ID lowercase and prepend 'ml' (for 'molecule'). """
    return "ml{}".format(_id).lower()


def make_pdb_name(_id, io):
    """ Make a name for a temporary PDB file either read 'in' or put 'out' from CHARMM. """
    if io not in ["in", "out"]:
        raise ValueError("The IO parameter must be 'in' or 'out,' not '{}.'".format(io))
    # All of the chains need a separate PDB file, so put the chain name at the end.
    return "{}chain{}.pdb".format(io, _id).lower()
