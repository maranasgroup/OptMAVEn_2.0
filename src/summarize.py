""" Create a table of summary data for experiments.
usage:
python summary.py summary_output_file.csv
"""

import csv
import os
import sys

import standards

output = sys.argv[1]
if os.path.exists(output):
    raise IOError("File exists: {}".format(output))

# fields in the output CSV file
fields = ["Experiment", "Chains", "Epitope Residues"]
# map keys in Summary.txt to fields in the CSV file
keys = {"Experiment name": "Experiment",
        "Antigen chains": "Chains",
        "Epitope residues": "Epitope Residues"}

experiments = list()
# Loop through all the experiments.
for experiment in os.listdir(standards.ExperimentsDirectory):
    exp_dir = os.path.join(standards.ExperimentsDirectory, experiment)
    summary_file = os.path.join(exp_dir, "Summary.txt")
    if not os.path.isfile(summary_file):
        # Skip the experiment if there is no Summary.txt file.
        continue
    with open(summary_file) as f:
        values = dict()
        for line in f:
            if ":" not in line:
                 # A ':' indicates a field; skip lines without a ':'
                 continue
            # Split the key from the value at the location of the first ':'
            key, value = line.split(":", 1)[0].strip(), "".join(line.split(":", 1)[1:]).strip()
            # Get the value of the key if the key is one of the keys to get.
            field = keys.get(key, None)
            if field == "Epitope Residues":
                # Add a space between epitope residues to improve readability.
                value = ", ".join([x.strip() for x in value.split(",")])
            if field is not None:
                values[field] = value
        # Add the information to the list of experiment information.
        experiments.append(values)

# Write the summary information.
with open(output, "w") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for experiment in sorted(experiments, key=lambda x: x["Experiment"]):
        writer.writerow(experiment)
