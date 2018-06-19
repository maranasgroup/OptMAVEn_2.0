""" The embedder module provides functions for embedding the MAPs parts in Euclidean space. """

from collections import OrderedDict
import itertools
import os
import pickle
import re
import shutil
import sys
import time

from Bio.SubsMat.MatrixInfo import blosum62 as b62
from Bio.SeqUtils import seq1
from Bio.Cluster import distancematrix
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as pl
import numpy as np

import maps
import standards


DGSolCommand = "dgsol"

gap = "_"

# Amino acid alphabet
aa = "ACDEFGHIKLMNPQRSTVWY" + gap

b62.update({(b, a): v for (a, b), v in b62.items()})  # make BLOSUM62 symmetric


def rank(array):
    order = array.argsort()
    ranks = np.empty(len(array), dtype=int)
    ranks[order] = np.arange(len(array))
    return ranks


def triangle_inequality(matrix):
    # NOTE: only works on square matrices
    if isinstance(matrix, dict):
        N = set()
        for (k1, k2) in matrix.keys():
            N.add(k1)
            N.add(k2)
        N = sorted(N)
    else:
        matrix = np.array(matrix)
        N = range(len(matrix))
    for i,j,k in itertools.product(N, repeat=3):
        if matrix[i,j] + matrix[j,k] < matrix[i,k]:
            return False
    return True


def get_residues(pdb_file):
    return OrderedDict([(int(line[22: 26]), seq1(line[16: 20].strip())) for line in open(pdb_file) if line.startswith("ATOM")])


def similarity(seq1, seq2, penalty_open, penalty_ext):
    resnums = sorted(set(seq1).union(set(seq2)))
    resis1 = "".join([seq1.get(i, gap) for i in resnums])
    resis2 = "".join([seq2.get(i, gap) for i in resnums])
    sim = 0
    open_gap = False
    for r1, r2 in zip(resis1, resis2):
        if r1 == gap or r2 == gap:
            if open_gap:
                sim -= penalty_ext
            else:
                sim -= penalty_open
                open_gap = True
        else:
            sim += b62[r1, r2]
            open_gap = False
    return sim


def distance_and_similarity(seq1, seq2, penalty_open, penalty_ext):
    s11 = similarity(seq1, seq1, penalty_open, penalty_ext)
    s12 = similarity(seq1, seq2, penalty_open, penalty_ext)
    s22 = similarity(seq2, seq2, penalty_open, penalty_ext)
    dist = max(s11 - s12, s22 - s12)
    return dist, s12

clustering_dir = os.path.join(standards.MapsDirectory, "clustering")
if not os.path.isdir(clustering_dir):
    os.mkdir(clustering_dir)

for gap_penalty in standards.MapsGaps:
    print("GAP PENALTY: {}".format(gap_penalty))
    gap_dir = os.path.join(clustering_dir, str(gap_penalty))
    if not os.path.isdir(gap_dir):
        os.mkdir(gap_dir)
    best_dir = os.path.join(gap_dir, "best_results")
    if not os.path.isdir(best_dir):
        os.mkdir(best_dir)

    coord_lines = list()  # store the best coordinates for all parts
    # Loop through each family (HV, HCDR3, HJ, LV, etc.) of MAPs parts.
    for chain, region in itertools.product(standards.MapsChains, standards.MapsRegions):
        category = chain + region
        print("CATEGORY: " + category)
        cat_dir = os.path.join(gap_dir, category)
        if not os.path.isdir(cat_dir):
            os.mkdir(cat_dir)
        # Dictionaries to store the sizes of the compressed contact maps.
        part_distances = OrderedDict()
    
        # Store part data.
        part_names = list()
        part_residues = dict()
        part_self_scores = dict()
        
        print("Listing MAPs parts")
        # For each part in the family:
        for part, part_file in sorted(maps.parts.items(), key=lambda x: int(maps.get_part_number(x[0]))):
            if part.startswith(category):
                part_names.append(part)
                # Read the residues.
                part_residues[part] = get_residues(part_file)
        N = len(part_names)

        distance_file = os.path.join(cat_dir, "{}{}_distances.csv".format(chain, region))
        if not os.path.isfile(distance_file) or True:#FIXME
            print("Computing MAPs alignment distances")
            Smat = np.zeros((N, N))
            # For each pair of parts:
            for i, part1 in enumerate(part_names):
                for j, part2 in enumerate(part_names[i:]):
                    # Compute the distance.
                    dist, sim = distance_and_similarity(part_residues[part1], part_residues[part2], gap_penalty, gap_penalty)
                    if part1 not in part_distances:
                        part_distances[part1] = OrderedDict()
                    part_distances[part1][part2] = dist
                    if part2 not in part_distances:
                        part_distances[part2] = OrderedDict()
                    part_distances[part2][part1] = dist
                    Smat[i, i + j] = sim
                    Smat[i + j, i] = sim
            """
            print("Validating MAPs alignment distances")
            # Ensure s criteria are met.
            for x, y in itertools.product(range(N), repeat=2):
                if Smat[x, x] >= Smat[x, y]:
                    pass
                else:
                    raise ValueError("Violated s({x}, {x}) >= s({x}, {y})".format(x=part_names[x], y=part_names[y]))
            print("Passed self-similarity")
            
            for x, y in itertools.product(range(N), repeat=2):
                if Smat[x, y] == Smat[y, x]:
                    pass
                else:
                    raise ValueError("Violated s({x}, {y}) = s({y}, {x})".format(x=part_names[x], y=part_names[y]))
            print("Passed symmetry")
        
            for x, y in itertools.product(range(N), repeat=2):
                if (Smat[x, y] == Smat[x, x] and Smat[y, x] == Smat[y, y]):
                    if x == y:
                        pass
                    else:
                        print([(part_residues[part_names[x]][ri], part_residues[part_names[y]][ri]) for ri in set(part_residues[part_names[x]]).union(part_residues[part_names[y]]) if part_residues[part_names[x]][ri] != part_residues[part_names[y]][ri]])
                        raise ValueError("Violated s({x}, {x}) = s({x}, {y}) = s({y}, {y}) ==> {x} = {y}".format(x=part_names[x], y=part_names[y]))
            print("Passed identity")
            
            for x, y, z in itertools.product(range(N), repeat=3):
                if Smat[x, y] + Smat[y, z] <= Smat[x, z] + Smat[y, y]:
                    pass
                else:
                    print(part_residues[part_names[x]], part_residues[part_names[y]])
                    raise ValueError("Violated s({x}, {y}) + s({y}, {z}) <= s({x}, {z}) + s({y}, {y})".format(x=part_names[x], y=part_names[y], z=part_names[z]))
            print("Passed similarity triangle")
            
            if not triangle_inequality([[part_distances[i][j] for j in part_names] for i in part_names]):
                raise ValueError("Violated triangle inequality.")
            print("Passed distance triangle")
            Dmat = np.array([[part_distances[i][j] for j in part_names] for i in part_names])
            N = Dmat.shape[0]
            e = np.ones((N, 1), dtype=np.float32)
            s = e / N
            f = (np.eye(N) - np.dot(e, s.T))
            F = np.dot(f, np.dot(Dmat, f))
            u, v = np.linalg.eigh(F)
            # Count the #s of eigenvalues >0 ("Pos"), <0 ("Neg"), and ==0 ("Zero").
            print("Pos", len([ui for ui in u if abs(ui) > 0.001 and ui > 0]), "Zero", len([ui for ui in u if abs(ui) <= 0.001]), "Neg", len([ui for ui in u if abs(ui) > 0.001 and ui < 0]))
            """
            
            # Write the distances to a CSV file.
            with open(distance_file, "w") as f:
                f.write(",".join([""] + part_names))
                for part1, parts in part_distances.iteritems():
                    f.write("\n" + ",".join([part1] + map(str, parts.values())))

        # Load the distances from the CSV file.   
        with open(distance_file) as f:
            parts = f.readline().strip().split(",")[1:]
            part_distances = OrderedDict([(line.split(",")[0], OrderedDict(zip(parts, map(float, line.split(",")[1:])))) for line in f])
                
        tol = 0.2
        incr = 0.05
        window = 0.0
        windows = list()
        spearmans = list()
        max_errors = list()

        report_file = os.path.join(cat_dir, "{}_report.txt".format(category))
        open(report_file, "w")

        # Test windows from 0 to 0.5 in increments of 0.05.
        while window <= 0.5:
            print("Window: {}".format(window))
            # Make the distance data file for dgsol.
            dgdata = os.path.join(cat_dir, "{}{}_win{}_dists.data".format(chain, region, window))
            open(dgdata, "w").write("\n".join(["{:>10} {:>9} {:>19} {:>19}".format(p1 + 1, p2 + 1, part_distances[part_names[p1]][part_names[p2]] * (1 - window), part_distances[part_names[p1]][part_names[p2]] * (1 + window)) for p1, p2 in itertools.combinations(range(len(part_names)), 2)]))   
            # Run dgsol.
            dgsol = os.path.join(cat_dir, "{}{}_win{}_coords.sol".format(chain, region, window))
            dgsum = os.path.join(cat_dir, "{}{}_win{}_coords.sum".format(chain, region, window))
            if not os.path.exists(dgsol) or not os.path.exists(dgsum) or len(open(dgsol).readlines()) != N + 3:
                print("Embedding MAPs alignment distances (window = {})".format(window))
                os.system("{} {} {} {}".format(DGSolCommand, dgdata, dgsol, dgsum))
            # Read the summary file.
            lines = open(dgsum).read().strip().split("\n")
            if len(lines) != 5 or len(lines[4].split()) != 6:
                # If there are not 5 lines in the file, something went wrong. Increase window.
                error = None
            else:
                error = float(lines[4].split()[5])
            
            # Read the coordinates.
            print("Saving MAPs coordinates (window = {})".format(window))
            coords = np.array([[float(x) for x in line.split()] for line in open(dgsol).readlines()[: N]])
            open(os.path.join(cat_dir, "{}{}_win{}_coords.csv".format(chain, region, window)), "w").write("\n".join(["{},{}".format(part, ",".join(map(str, coord))) for part, coord in zip(part_names, coords)]))
            # Compare the optimal distances with the target distances.
            D_target = np.array([[part_distances[i][j] for j in part_names] for i in part_names])
            print("Computing pairwise coordinate distances (window = {})".format(window))
            D_optimized = np.array([[np.linalg.norm(coords[i] - coords[j]) for j in range(N)] for i in range(N)])
            open(os.path.join(cat_dir, "{}{}_win{}_dists.csv".format(chain, region, window)), "w").write("\n".join([",".join(map(str, row)) for row in D_optimized]))
            dist_error = np.max(np.abs(D_target - D_optimized))
            print("Flattening matrices (window = {})".format(window))
            D_target_flat = np.hstack([D_target[i, i:] for i in range(N)]).flatten()
            D_opt_flat = np.hstack([D_optimized[i, i:] for i in range(N)]).flatten()
            print("Calculating spearman (window = {})".format(window))
            spearman = np.corrcoef(rank(D_target_flat), rank(D_opt_flat))[0, 1]
            xmax, ymax = max(D_target_flat), max(D_opt_flat)
            figure = os.path.join(cat_dir, "{}{}_win{}_dists.png".format(chain, region, window))
            print("Making figure (window = {})".format(window))
            pl.plot(D_target_flat, D_opt_flat, "o", c="gray", alpha=0.1 + 0.9 * 2.71**(-N / 50.0), mew=0.0)
            pl.xlabel("Pairwise alignment distance")
            pl.ylabel("Euclidean distance of coordinates")
            pl.title(chain + region)
            pl.text(xmax * 0.15, ymax * 0.75, "$\\rho = {:.5f}$".format(spearman), fontsize=18)
            pl.savefig(figure)
            pl.close()
            open(report_file, "a").write("{:>10} {:>9.5f} {:>9.5f}\n".format(window, dist_error, spearman))
            max_errors.append(dist_error)
            spearmans.append(spearman)
            windows.append(window)        
            window += incr
        print("Finding optimal window")
        order = np.argsort(spearmans)
        spearmans = np.array(spearmans)[order]
        max_errors = np.array(max_errors)[order].compress(spearmans == spearmans[-1])
        windows = np.array(windows)[order].compress(spearmans == spearmans[-1])
        window = windows[np.argmin(max_errors)]
        for f in ["coords.csv", "dists.png"]:
            fname = "{}{}_win{}_{}".format(chain, region, window, f)
            shutil.copyfile(os.path.join(cat_dir, fname), os.path.join(best_dir, fname))
        coord_lines.extend([line.strip() for line in open(os.path.join(best_dir, "{}{}_win{}_coords.csv".format(chain, region, window))).readlines()])
    # Copy the best coordinates into the MAPs directory.
    coord_file = os.path.join(standards.MapsDirectory, "coords_{}.csv".format(gap_penalty))
    if "replace" in sys.argv or not os.path.isfile(coord_file):
        with open(coord_file, "w") as f:
            f.write("\n".join(coord_lines))
