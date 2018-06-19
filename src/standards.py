""" This module stores information that is shared between the modules. """

from collections import OrderedDict
import itertools
import os
import random
import shutil
import string
import sys

import Bio
import numpy as np


# Environment
PythonCommand = ""
CharmmCommand = ""
CplexDirectory = ""
VmdCommand = ""
NamdCommand = ""
PbsQueue = ""

# Versioning
OptmavenVersion = "2.0"
OptmavenName = "OptMAVEn-{}".format(OptmavenVersion)

SupportedPythonVersion = "2.7"
if not sys.version.startswith(SupportedPythonVersion):
    raise NotImplementedError("{} requires Python {}, got {}".format(OptmavenName, SupportedPythonVersion, sys.version))
SupportedBiopythonVersion = "1.7"
if not Bio.__version__.startswith(SupportedBiopythonVersion):
    raise NotImplementedError("{} requires Biopython {}, got {}".format(OptmavenName, SupportedBiopythonVersion, Bio.__version__))

SourceDirectory = os.path.dirname(os.path.realpath(__file__))
OptmavenDirectory = os.path.dirname(SourceDirectory)
ExperimentsDirectory = os.path.join(OptmavenDirectory, "experiments")
DataDirectory = os.path.join(OptmavenDirectory, "data")
PDBDirectory = os.path.join(DataDirectory, "pdbs")
InputsDirectory = os.path.join(DataDirectory, "input_files")
AntibodiesDirectory = os.path.join(DataDirectory, "antibodies")
ScaffoldsDirectory = os.path.join(AntibodiesDirectory, "scaffolds")
MapsDirectory = os.path.join(AntibodiesDirectory, "maps")

# Check to make sure all of the directories and files are present. If not, the installation is malformed.
if os.path.basename(OptmavenDirectory) != OptmavenName:
    raise OSError("This installation of Optmaven is located in {} but needs to be located in {}.".format(OptmavenDirectory, OptmavenName))

# Make the experiments directory.
if not os.path.isdir(ExperimentsDirectory):
    os.mkdir(ExperimentsDirectory)

# Maximum line length in console.
ConsoleWidth = 80
# Maximum number of list items to print.
MaxListDisplay = 20
SelectAll = "all"
SelectNone = "none"
EscapeCharacter = "\\"

# Datetime settings.
DatetimeFormat = "%Y %b %d %H:%M:%S"

# Path configuration.
AllowedPathCharacters = string.letters + string.digits + "_."

# Experimental configuration.
DefaultNumberOfDesigns = 5000

# Argument standards.
HetAtmExclude = "yes"
HetAtmInclude = "no"
HetAtmAsk = "ask"
HetAtmOptions = [HetAtmAsk, HetAtmInclude, HetAtmExclude]

# Default PBS settings.
DefaultWalltime = 86399  # 23 hours, 59 minutes, 59 seconds
DefaultBatchSize = 1
PbsQsub = "qsub"
PbsArrayId = "$PBS_ARRAYID"
PbsJobFilePrefix = "job-"
UnixTimeCodes = OrderedDict([("Real", "%e"), ("User", "%U"), ("System", "%S")])
PbsTimeFormat = "'{}'".format(r"\n".join(UnixTimeCodes.values()))
TimeCommand = "/usr/bin/time"

# Benchmarking settings.
BenchmarkingFields = ["Type", "Purpose", "Detail"] + UnixTimeCodes.keys() + ["Drive Usage", "Time Stamp"]

# Optmaven grid settings.
DefaultOptmavenGrid_x = np.linspace(-10, 5, 7)
DefaultOptmavenGrid_y = np.linspace(-5, 10, 7)
DefaultOptmavenGrid_z = np.linspace(3.75, 16.25, 11)
DefaultOptmavenGrid_zAngle = np.linspace(0, 300, 6)

# Input files.
DefaultTopologyFile = os.path.join(InputsDirectory, "top_all27_prot_na.rtf")
DefaultParameterFile = os.path.join(InputsDirectory, "par_all27_prot_na.prm")
DefaultSolvationFile = os.path.join(InputsDirectory, "solvation.dat")

# CHARMM configuration.
CharmmSolvationTerm = "gbener"
DefaultCharmmEnergyTerms = ['angl', 'bond', 'dihe', 'elec', 'impr', 'urey', 'vdw', CharmmSolvationTerm]
DefaultCharmmIterations = 5000

# VMD configuration.
VmdDisp = "-dispdev"
VmdNoDisp = "none"
VmdExec = "-e"
VmdMolecules = "-m"
VmdFrames = "-f"
VmdArgs = "-args"
VmdFunctions = os.path.join(SourceDirectory, "vmd_functions.tcl")

# Optmaven files.
ScaffoldChains = ["H", "K"]
ScaffoldAntibodies = {chain: os.path.join(ScaffoldsDirectory, "Molecule{}.pdb".format(chain)) for chain in ScaffoldChains}

# Antigen positioning parameters.
DefaultClashCutoff = 1.25  # Angstroms

# MAPs database standards.
MapsGaps = [4, 6, 8, 10, 12]
MapsHeavyChains = ["H"]
MapsLightChains = ["K", "L"]
MapsNamesakeHeavy = "H"
MapsNamesakeLight = "L"
MapsChains = MapsHeavyChains + MapsLightChains
MapsNamesakeChains = [MapsNamesakeHeavy, MapsNamesakeLight]
MapsRegions = ["V", "CDR3", "J"]
MapsCdrs, MapsNamesakeCdrs = [["{}{}".format(chain, cdr) for chain, cdr in itertools.product(chains, MapsRegions)] for chains in [MapsChains, MapsNamesakeChains]]
MapsIntegerCutsFile = os.path.join(MapsDirectory, "MAPs_Integer_Cuts.txt")

# Antibody coordinate standards.
MapsCoordSep = "."
MapsCoordDimension = 3

def make_maps_coord(cdr, dimension):
    """ Convert a CDR and a dimension to a string.
    e.g. 'LV' and 2 becomes 'LV.2'
    Each MAPs part is a 3-dimensional vector; the dimensions for e.g. LV are 'LV.1', 'LV.2', 'LV.3'
    """
    return "{}{}{}".format(cdr, MapsCoordSep, dimension)

def to_maps_coord(item):
    """ Convert a string back into a CDR and dimension. """
    if MapsCoordSep in item:
        split = item.split(MapsCoordSep)
        if len(split) == 2:
            cdr, dimension = split
            if cdr in MapsCdrs and str(dimension) in map(str, range(MapsCoordDimension)):
                return cdr, int(dimension)
    raise TypeError("Cannot convert to MAPs coord: {}".format(item))

AngleLabel = "Angle"
xLabel = "x"
yLabel = "y"
zLabel = "z"
zAngleLabel = zLabel + AngleLabel
SinLabel = "Sin"
CosLabel = "Cos"
PositionOrder = [zAngleLabel, xLabel, yLabel, zLabel]
AtomCoordOrder = [xLabel, yLabel, zLabel]
PositionCoordOrder = list(itertools.chain(*[["{}{}".format(item, label) for label in [SinLabel, CosLabel]] if AngleLabel in item else [item] for item in PositionOrder]))
MapsCoordOrder = [make_maps_coord(cdr, dim) for cdr, dim in itertools.product(MapsNamesakeCdrs, range(MapsCoordDimension))]
EnergyLabel = "energy"
CoordOrder = PositionCoordOrder + MapsCoordOrder
DefaultGapPenalty = 8

# Kmeans standards.
DefaultKmeansMaxIterations = 1000
DefaultKmeansTolerance = 0.01
DefaultKmeansOptKThreshold = 0.25

# Rotation matrix routines.

_ROTATION_MATRIX_DIMENSION = len(AtomCoordOrder)
dim = _ROTATION_MATRIX_DIMENSION
xAxis, yAxis, zAxis = np.eye(dim)

def degrees_to_radians(degrees):
    return degrees * np.pi / 180.0


def radians_to_degrees(radians):
    return radians * 180.0 / np.pi


def rotate_vi_to_vf(vi, vf):
    """ Compute a rotation matrix to rotate the 3D coordinate vi to vf. """
    if len(vi) != dim or len(vf) != dim:
        raise ValueError("The rotation matrix function requires 3D coordinates.")
    # Normalize both vectors.
    vi, vf = np.array(vi), np.array(vf)
    vin = np.linalg.norm(vi)
    vfn = np.linalg.norm(vf)
    if np.isclose(vin, 0) or np.isclose(vfn, 0):
        # Abort rotation if either vector is almost zero.
        return np.eye(dim)
    vi /= vin
    vf /= vfn
    # Find the axis of rotation, which is perpendicular to both vectors.
    ax = np.cross(vf, vi)
    # Calculate the sine and cosine of the angle of rotation.	
    sin = np.linalg.norm(ax)
    cos = np.dot(vf, vi)
    return rotate_axis_angle(ax, sin=sin, cos=cos)


def rotate_axis_angle(axis, angle=None, degrees=True, sin=None, cos=None):
    """ Create a rotation matrix to rotate by an angle around an axis. """
    if len(axis) != dim:
        raise ValueError("The rotation matrix function requires 3D coordinates.")
    if angle is not None:
        if degrees:
            angle = degrees_to_radians(angle)
    if sin is None:
        sin = np.sin(angle)
    if cos is None:
        cos = np.cos(angle)
    axn = np.linalg.norm(axis)
    if np.isclose(axn, 0):
        # If the rotation axis is almost zero, then the rotation angle is either almost 0 or 180.
        # If the angle is 0, then cos = 1. Return identity matrix.
        # If the angle is 180, then cos = -1. Return the negative identity matrix.
        return cos * np.eye(dim)
    # Normalize the axis.
    ax = axis / axn
    # Compute the cross product matrix of the axis.
    cpm = np.array([
    	[     0, -ax[2],  ax[1]],
    	[ ax[2],      0, -ax[0]],
    	[-ax[1],  ax[0],      0]
    ])
    # Compute the tensor product matrix of the axis.
    tpm = ax.reshape(1, dim) * ax.reshape(dim, 1)
    # Return the rotation matrix.
    return cos * np.eye(dim) + sin * cpm + (1 - cos) * tpm


def random_string(length, alphabet=None):
    if alphabet is None:
        alphabet = string.letters
    return "".join([random.choice(alphabet) for i in range(length)])


def is_path_component(file_name):
    """ Determine if a file name is a single file without whitespace. """
    return isinstance(file_name, str) and len([char for char in file_name if char not in AllowedPathCharacters]) == 0


def is_subdirectory(directory, parent_directory):
    directory = os.path.realpath(directory)
    parent_directory = os.path.realpath(parent_directory)
    parent_path_old = directory
    parent_path_new, child = os.path.split(parent_path_old)
    while parent_path_old != parent_path_new:
        if parent_path_new == parent_directory:
            return True
        parent_path_old = parent_path_new
        parent_path_new, child = os.path.split(parent_path_old)
    return False


# Safe removal of directory trees.
def safe_rmtree(directory):
    if is_subdirectory(directory, ExperimentsDirectory):
        shutil.rmtree(directory)
    else:
        raise OSError("Directory trees may only be removed if they are subdirectories of {}".format(ExperimentsDirectory))


# Safe removal of entire experiments.
def safe_rm_experiment(experiment):
    full_experiment = os.path.join(ExperimentsDirectory, experiment)
    if is_path_component(experiment) and experiment in os.listdir(ExperimentsDirectory) and os.path.isdir(full_experiment):
        if raw_input("Are you SURE you want to remove the Experiment {}? You CANNOT undo this action. ".format(experiment)).lower().startswith("y"):
            shutil.rmtree(full_experiment)
            if os.path.isdir(full_experiment):
                print("Optmaven could not completely remove {}.".format(experiment))
            else:
                print("The Experiment {} has been removed.".format(experiment))
        else:
            print("The Experiment {} was not removed.".format(experiment))
    else:
        raise OSError("Cannot remove Experiment: {} does not exist.".format(full_experiment))
