""" Print the current status of each experiment on the screen.

usage:
python check_status [prefix]

If prefix is given, then only experiments that begin with the prefix will be listed.
"""


import cPickle as pkl
import os
import sys

import console
import standards


console.clear()
console.disp("{:<39} {:<39}".format("EXPERIMENT", "STATUS"))
console.disp("-" * 80)
for experiment in os.listdir(standards.ExperimentsDirectory):
    if len(sys.argv) > 1:
        if not any([x in experiment for x in sys.argv[1:]]):
            continue
    directory = os.path.join(standards.ExperimentsDirectory, experiment)
    errors = os.path.join(directory, "errors.txt")
    pickle = os.path.join(directory, ".temp", "{}.pickle".format(experiment))
    summary = os.path.join(directory, "Summary.txt")
    results = os.path.join(directory, "Results.csv")
    if os.path.isfile(errors):
        status = "ERROR: see errors.txt"
    else:
        if all([os.path.isfile(x) for x in [summary, results]]):
            status = "Completed"
        else:
            try:
                with open(pickle) as f:
                    exp = pkl.load(f)
            except IOError:
                status = "ERROR: missing pickle, results, summary, and error files"
            else:
                task, status = exp.get_task(exp.status)
    console.disp("{:<39} {}".format(experiment, status))
