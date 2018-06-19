""" Make a pickle file that indexes all of the MAPs parts as a dict of
{category_1: {part_1: path_1, part_2, path_2, ... },
 category_2: {part_1: path_1, part_2, path_2, ... }, ... } """

from collections import OrderedDict
import cPickle as pkl
import itertools
import os
import re

domains = ["H", "K", "L"]
regions = ["V", "CDR3", "J"]

MAPs_directory = os.path.dirname(os.path.realpath(__file__))
parts = dict()
for domain, region in itertools.product(domains, regions):
    category = domain + region
    parts[category] = OrderedDict([(os.path.splitext(pdb)[0], os.path.join(
            MAPs_directory, category, pdb)) for pdb in os.listdir(os.path.join(
            MAPs_directory, category)) if re.match(category + "_[0-9]+.pdb",
            pdb)])
pkl.dump(parts, open(os.path.join(MAPs_directory, "MAPs_Index.pickle"), "w"))
