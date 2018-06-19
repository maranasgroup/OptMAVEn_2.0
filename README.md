Optmaven-2.0
16 June 2018
Written by Matthew F. Allan and Ratul Chowdhury with contributions from Tong Li and Robert J. Pantazes
Costas Maranas Laboratory, Department of Chemical Engineering
The Pennsylvania State University
University Park, PA 16802


DEPENDENCIES

Python 2.7: https://www.python.org/
NumPy 1.13+: http://www.numpy.org/
SciPy 0.19+: https://www.scipy.org/
Biopython 1.7+: https://biopython.org/
VMD 1.9+: http://www.ks.uiuc.edu/Research/vmd/
NAMD 2.12+: http://www.ks.uiuc.edu/Research/namd/
CHARMM 34+: https://www.charmm.org/
CPLEX - Python API: https://www.ibm.com/support/knowledgecenter/SSSA5P_12.7.1/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/setup_overview.html
DGSOL (only if running embedder module): http://www.mcs.anl.gov/~more/dgsol/


INSTALLATION

OptMAVEn-2.0 needs no installation. Merely unzip and/or untar the compressed
directory and move it to the location of choice. However, you must edit the
file src/standards.py by specifying the following environment variables:

PythonCommand: the path to the Python executable
CharmmCommand: the path to the CHARMM executable
CplexDirectory: the directory in which the CPLEX modules are located
VmdCommand: the path to the VMD executable
NamdCommand: the path to the NAMD executable
PbsQueue: the name of the PBS queue to which to submit jobs

For example, change the line
PythonCommand = ""
to
PythonCommand = "/usr/bin/python"


DIRECTORY STRUCTURE

The OptMAVEN-2.0 directory contains the following subdirectories:

src: source modules in Python and TCL
data: data files needed for OptMAVEn-2.0 to run. Subdirectories:
    antibodies: MAPs database and framework antibody structures
    input_files: CHARMM topology, parameter, and solvation files
    pdbs: structures of the input antigens
experiments: contains all experiment directories

The OptMAVEn-2.0 directory also contains the following executables:

OptMAVEN-2.0: start an OptMAVEn-2.0 experiment
find_contacts: find the residues that make contacts between two molecules
interaction_energy: compute the interaction energy between two molecules
check_status: check the status of the experiments
remove_experiment: remove an experiment permanently


RUNNING OPTMAVEN-2.0

Navigate to the main directory of OptMAVEn-2.0 and type ./OptMAVEn-2.0 to start
OptMAVEn-2.0. You will be prompted for the following information:

Name of the experiment
Whether or not to use benchmarking
The gap penalty to use during clustering
Whether or not to customize the following settings:
    Antigen positioning grid
    Antigen positioning clash cutoff
    CHARMM topology files
    CHARMM parameter files
    CHARMM solvation files
    CHARMM energy terms
    CHARMM iteration limit
    Queue walltime
    Queue batch size
File or PDB id of the antigen
Whether or not to exclude heteroatoms from each chain
The chain(s) that constitute(s) the antigen
The residue(s) that constitute(s) the epitope to target

The following command-line arguments may be specified for OptMAVEn-2.0:
--keeptemp: keep temporary files. Useful for debugging but requires more space
--exclude_hetero: 'ask' (default), 'yes', or 'no'
    'ask': ask if heteroatoms should be excluded
    'yes': automatically exclude all heteroatoms
    'no': do not exclude any heteroatoms


CHECKING THE STATUS

The status may be checked by running ./check_status from the main directory.
Supplying a command-line argument returns only the experiments whose names
contain the argument. For example:
./check_status exp
would report the status of experiments named 'exp', 'exp2', and 'my_experiment'
but not 'myExperiment'.


VIEWING THE RESULTS

The results of an experiment 'EXP_NAME' are saved in the file
experiments/EXP_NAME/Results.csv. Additionally, a summary of the inputs is given
in experiments/EXP_NAME/Summary.txt.


TROUBLESHOOTING

If an experiment 'EXP_NAME' fails, an error message will be written to the file
experiments/EXP_NAME/errors.txt. Please review the error message to diagnose
the problem. In the case that no error message is generated, the problem likely
arose from the PBS queueing system. Please ensure that the queue is set up
correctly and verify all environment settings (e.g. the path to the Python
executable defined in src/standards.py).