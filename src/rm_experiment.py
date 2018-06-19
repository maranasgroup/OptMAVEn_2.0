""" Safely remove an experiment's directory.
usage:
python rm_experiment.py experiment_name
"""

import sys

import standards

experiment = sys.argv[1]
standards.safe_rm_experiment(experiment)
