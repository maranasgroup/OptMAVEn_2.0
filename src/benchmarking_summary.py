""" Create a summary of the benchmarking results of an experiment.

usage:
python benchmarking_summary.py name_of_summary_report.csv

The file name_of_summary_report.csv will be created.
It will have the fields given in FIELDS.
Experiment: experiment name
Real, User, System: sums of the real, user, and system times
Drive Usage: maximum drive usage measured
Chains, Residues, Atoms: the numbers of chains, residues, and atoms in the antigen
Positions, Selected Designs: the number of non-clashing positions, and the number of designs output at the end (Positions >= Selected Designs)
Min Unrelaxed Energy, Min Relaxed Energy: the antibody-antigen interaction energies (kcal/mol) of the antibody-antigen complex before and after relaxation
"""

import csv
import os
import sys

from Bio.PDB import PDBParser
import numpy as np

import standards


KEY_TOTAL = "Total"
KEY_TYPE = "Type"
# The fields 
FIELDS = ["Experiment", "Real", "User", "System", "Drive Usage", "Chains", "Residues", "Atoms", "Positions", "Selected Designs", "Min Unrelaxed Energy", "Min Relaxed Energy"]

# Get the name of the output file.
try:
    output = sys.argv[1]
    if os.path.exists(output):
        raise IOError("File exists: {}".format(output))
except IndexError:
    raise IOError("You must specify an output file.")


# Get benchmarking results from each experiment.
results = list()
for experiment in os.listdir(standards.ExperimentsDirectory):
    print(experiment)
    edir = os.path.join(standards.ExperimentsDirectory, experiment)
    bmark_file = os.path.join(edir, "Benchmarking.csv")
    summary_file = os.path.join(edir, "Summary.txt")
    results_file = os.path.join(edir, "Results.csv")
    # Skip experiments without benchmarking information.
    if not all([os.path.isfile(f) for f in [bmark_file, summary_file, results_file]]):
        continue
    with open(bmark_file) as f:
        for row in csv.DictReader(f):
            if row[KEY_TYPE] == KEY_TOTAL:
                results.append({k: v for k, v in row.items() if k in FIELDS})
                break
        else:
            raise ValueError("Cannot find {} row in file: {}".format(KEY_TOTAL, bmark_file))
    with open(summary_file) as f:
        positions, selected = None, None
        for line in f:
            if line.startswith("Total positions:"):
                positions = int(line.split(":")[1])
            elif line.startswith("Selected designs:"):
                selected = int(line.split(":")[1])
        if any([x is None for x in (positions, selected)]):
            raise ValueError("Cannot find all information in file: {}".format(summary_file))
    with open(results_file) as f:
        rows = sorted([row for row in csv.DictReader(f)], key = lambda x: float(x["Relaxed energy (kcal/mol)"]))
        # Exclude rows that lie outside of the [-5, 5] range in x and y coordinates.
        rows = [row for row in rows if abs(float(row["x"])) <= 5 and abs(float(row["y"])) <= 5]
        try:
            min_relaxed_energy = float(rows[0]["Relaxed energy (kcal/mol)"])
            min_unrelaxed_energy = float(rows[0]["Unrelaxed energy (kcal/mol)"])
        except IndexError:
            min_unrelaxed_energy = np.nan
            min_relaxed_energy = np.nan
    # Count the number of atoms and residues in the antigen.
    ag_file = os.path.join(standards.ExperimentsDirectory, experiment, "structures", "antigen_relaxed.pdb")
    structure = PDBParser().get_structure("antigen", ag_file)
    n_chains, n_atoms, n_residues = 0, 0, 0
    for chain in structure.get_chains():
        n_chains += 1
        n_residues += len(list(chain.get_residues()))
        n_atoms += len(list(chain.get_atoms()))
    results[-1].update({"Experiment": experiment, "Chains": n_chains, "Residues": n_residues, "Atoms": n_atoms, "Positions": positions, "Selected Designs": selected, "Min Unrelaxed Energy": min_unrelaxed_energy, "Min Relaxed Energy": min_relaxed_energy})

# Output results as a CSV file.
with open(output, "w") as f:
    writer = csv.DictWriter(f, FIELDS)
    writer.writeheader()
    for experiment in sorted(results):
        writer.writerow(experiment)

