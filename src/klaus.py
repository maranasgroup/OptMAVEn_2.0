__doc__ = """ This module is the interface between IPRO and Klaus Schulten
programs VMD and NAMD and contains a built-in function for performing each task
 with VMD and NAMD. """

from collections import OrderedDict 
import itertools
import os
import shutil
import subprocess
import tempfile
import time

import numpy as np

import maps
import molecules
import charmm
import standards
import submitter

# Define how to run VMD.
vmd_command = standards.VmdCommand


class VmdProc(object):
    """ This class is the base class for managing all VMD procedures. """
    def __init__(self, experiment, directory=None):
        self.experiment = experiment
        self.directory = directory
        # Indicate that VMD has not yet been run.
        self.has_run = False

    def make_command(self, output=True):
        """ Create a VMD command. """
        # Start by listing all of the components common to all VMD commands:
        # VmdCommand: the path to the VMD executible
        # VmdDisp: a flag indicating which display mode to use
        # VmdNoDisp: an option for VmdDisp telling VMD to create no GUI
        # VmdExec: a flag telling VMD to execute the script that follows
        # script: the TCL script to run        
        command = [standards.VmdCommand, standards.VmdDisp, standards.VmdNoDisp, standards.VmdExec, os.path.join(standards.SourceDirectory, self.script)]
        # If VMD should load any molecules, load them by specifying their file paths after the VmdMolecules flag.
        try:
            command.append("{} {}".format(standards.VmdMolecules, " ".join(map(str, self.molecules))))
        except AttributeError:
            pass
        # If VMD should load any molecules as frames, load them by specifying their file paths after the VmdFrames flag.
        try:
            command.append("{} {}".format(standards.VmdFrames, " ".join(map(str, self.frames))))
        except AttributeError:
            pass
        # If VMD needs any additional arguments, specify them after the VmdArgs flag.
        try:
            command.append("{} {}".format(standards.VmdArgs, " ".join(map(str, self.args))))
        except AttributeError:
            pass
        # Optionally, redirect the output to vmd.out and vmd.err using the redirection symbols '>' and '2>'
        if output:
            outfile = os.path.join(self.directory, "vmd.out")
            errfile = os.path.join(self.directory, "vmd.err")
            command.extend([">", outfile, "2>", errfile])
        return command

    def vmd(self):
        """ Run VMD. """
        # List the command arguments.
        arguments = self.make_command()
        # Join them into a single command.
        command = " ".join(arguments)
        # Run the command.
        status = os.system(command)
        # Indicate if there was non-zero exit status.
        if status != 0:
            raise RuntimeError("Running VMD with the following command has failed:\n{}".format(command))

    def collect_garbage(self):
        """ Remove all temporary files for the process. """
        if not self.experiment.args.keeptemp:
            # Only remove the files if keeptemp has not been specified.
            try:
                os.remove(self.vmd_functions_file)
            except OSError:
                pass
            self.experiment.safe_rmtree(self.directory)

    def write_epitope_file(self):
        """ Write a file specifying the epitope residues. """
        self.epitope_file = os.path.join(self.directory, "epitope.txt")
        with open(self.epitope_file, "w") as f:
            f.write(atomselect_residues(self.experiment.epitope_residue_ids))

    def __enter__(self):
        # Create a directory for the process.
        if self.directory is None:
            self.directory = tempfile.mkdtemp(prefix="vmd_", dir=self.experiment.temp)
        else:
            os.mkdir(self.directory)
        # Make a link to the VMD functions file.
        self.vmd_functions_file = os.path.join(self.directory, os.path.basename(standards.VmdFunctions))
        os.symlink(standards.VmdFunctions, self.vmd_functions_file)
        # VMD must be run from the directory of this process.
        # So that the program can return to the original directory after VMD finishes, record the original directory.
        self.previous_directory = os.getcwd()
        # Change to the process directory.
        os.chdir(self.directory)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # Change back to the original directory.
        os.chdir(self.previous_directory)
        # Remove all temporary files.
        self.collect_garbage()


class PositionAntigen(VmdProc):
    """ Move an antigen to a new position. """
    def __init__(self, experiment, input_file, output_file, zAngle, x, y, z):
        VmdProc.__init__(self, experiment)
        # The input molecule, output file, antigen position, and experiment must all be specified.
        self.molecules = [input_file]
        self._args = [experiment.antigen_chain_ids[0], output_file, zAngle, x, y, z]

    def __enter__(self):
        VmdProc.__enter__(self)
        # The position is defined by the epitope center of geometry. Thus, the epitope residues must be known.
        self.write_epitope_file()
        # Add the path to the epitope file to the arguments.
        self.args = [self.epitope_file] + self._args
        # Run the script in VMD.
        self.script = "position_antigen.tcl"
        self.vmd()
        return self


class GridSearch(VmdProc):
    """ Generate non-clashing antigen positions using a grid search. """
    def __init__(self, experiment):
        VmdProc.__init__(self, experiment)

    def write_grid_file(self):
        """ Write a file specifying the levels of x, y, z, and zAngle to search. """
        self.grid_file = os.path.join(self.directory, "grid.dat")
        with open(self.grid_file, "w") as f:
            for dimension, levels in [
                (standards.xLabel, self.experiment.grid_x),
                (standards.yLabel, self.experiment.grid_y),
                (standards.zLabel, self.experiment.grid_z),
                (standards.zAngleLabel, self.experiment.grid_zAngle)
            ]:
                f.write("{}: {}\n".format(dimension, " ".join(map(str, levels))))

    def __enter__(self):
        VmdProc.__enter__(self)
        # First, write the grid file so that VMD knows which positions to search.
        self.write_grid_file()
        # The positions are defined by the epitope center of geometry.
        self.write_epitope_file()
        # Define the output file of non-clashing positions.
        self.positions_file = os.path.join(self.directory, "positions.dat")
        # The antigen will be tested for clashes with the two scaffold antibodies.
        self.molecules = [self.experiment.epitope_zmin_file, standards.ScaffoldAntibodies["H"], standards.ScaffoldAntibodies["K"]]
        # The arguments must go in this order.
        self.args = [self.grid_file, self.epitope_file, self.experiment.antigen_chain_ids[0], self.positions_file, self.experiment.clash_cutoff]
        # Run VMD.
        self.script = "grid_search.tcl"
        self.vmd()
        # Move the positions file from the temporary directory (which will be deleted when this process exits) to its proper place in the experiment directory.
        shutil.move(self.positions_file, self.experiment.positions_file)
        return self
        

class MergeAntigenMaps(VmdProc):
    """ Calculate the interaction energies between the antigen and a MAPs part in every position. """
    def __init__(self, experiment, maps_part, destination_directory):
        VmdProc.__init__(self, experiment)
        # The MAPs part must be specified.
        self.part = maps_part
        # The directory in which the file of energies will be located.
        self.destination_directory = destination_directory

    def __enter__(self):
        VmdProc.__enter__(self)
        # Merge the antigen and antibody into one molecule and create a PSF.
        # This step is necessary for NAMDEnergy, which needs a PSF and can only calculate the interaction energy between atoms in the same molecule.
        self.prefix = os.path.join(self.directory, "merged")
        pdb = "{}.pdb".format(self.prefix)
        psf = "{}.psf".format(self.prefix)
        # These arguments are required to create the PSF.
        self.args = [self.experiment.epitope_zmin_file, maps.parts[self.part], self.prefix]
        self.args.extend(self.experiment.topology_files)
        # Create the PSF and merge the molecules with VMD.
        self.script = "merge_antigen_part.tcl"
        self.vmd()
        # Move the PDB and PSF files of the merged molecule from the temporary directory (which will be deleted after this process exits) to their proper places. 
        shutil.move(pdb, self.destination_directory)
        shutil.move(psf, self.destination_directory)
        # Record the paths of the PDB and PSF files so that the function that created this process has access to these paths.
        self.pdb = os.path.join(self.destination_directory, os.path.basename(pdb))
        self.psf = os.path.join(self.destination_directory, os.path.basename(psf))
        return self


class MapsEnergies(VmdProc):
    """ Calculate the interaction energies between the antigen and a MAPs part in every position. """
    def __init__(self, experiment, maps_part, energy_finished):
        VmdProc.__init__(self, experiment)
        # The MAPs part must be specified.
        self.part = maps_part
        # The final location to which to move the file of energies must also be specified.
        self.energy_finished = energy_finished

    def __enter__(self):
        VmdProc.__enter__(self)
        # NAMDEnergy requires a PDB and PSF file of the merged antigen and antibody.
        with MergeAntigenMaps(self.experiment, self.part, self.directory) as merge:
            self.pdb = merge.pdb
            self.psf = merge.psf
        # Name a temporary file in which to record the energies as they are calculated.
        energy_temp = os.path.join(self.directory, "energies_temp.dat")
        # Write the file specifying the epitope so that the antigen can be positioned properly.
        self.write_epitope_file()
        # These arguments are required.
        self.args = [self.psf, self.pdb, self.experiment.positions_file, self.epitope_file, self.experiment.antigen_chain_ids[0], energy_temp, standards.NamdCommand]
        self.args.extend(self.experiment.parameter_files)
        # Run VMD.
        self.script = "interaction_energies.tcl"
        self.vmd()
        # Move the temporary energy file to its permanent location.
        try:
            os.makedirs(os.path.dirname(self.energy_finished))
        except OSError:
            pass
        shutil.move(energy_temp, self.energy_finished)
        return self


class CreateAntibody(VmdProc):
    """ Merge six MAPs parts into an antibody. """
    def __init__(self, experiment, maps_parts, output_file):
        VmdProc.__init__(self, experiment)
        # The names of the six MAPs parts are required.
        self.maps_parts = dict()
        # The name of the output file of merged parts is also required.
        self.args = [output_file]
        # Ensure that exactly one of each V, CDR3, and J region is present for each of the light and heavy chains.
        for part in maps_parts:
            cdr, number = maps.split(part)
            self.maps_parts[maps.translate_chain_namesake(cdr)] = maps.parts[part]
        if sorted(self.maps_parts) != sorted(standards.MapsNamesakeCdrs):
            raise ValueError("An antibody needs one of each CDR, not: {}".format(", ".join(self.maps_parts)))
    
    def __enter__(self):
        VmdProc.__enter__(self)
        # Make sure the MAPs parts are loaded as molecules.
        self.molecules = [self.maps_parts[cdr] for cdr in standards.MapsNamesakeCdrs]
        # Run VMD.
        self.script = "create_antibody.tcl"
        self.vmd()
        return self


def atomselect_residues(residue_ids):
    """ Create a residue selection string in VMD atomselect syntax. """
    return " or ".join(["(chain {} and resid {})".format(c_id, "{}{}{}".format(*r_id)) for c_id, r_ids in residue_ids.items() for r_id in r_ids])
