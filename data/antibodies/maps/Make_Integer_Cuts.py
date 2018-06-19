""" Make an integet cut file of all clasing MAPs part pairs. """

from collections import OrderedDict
import cPickle as pkl
import itertools
import os
import re
import sys

MAPs_dir = os.path.dirname(os.path.realpath(__file__))
database_dir = os.path.dirname(MAPs_dir)
OptMAVEn_dir = os.path.join(database_dir, "OptMAVEn")
main_dir = os.path.dirname(database_dir)
mods_dir = os.path.join(main_dir, "modules")
sys.path.append(mods_dir)
import KS

if len(sys.argv) == 1:
    cut_file = os.path.join(OptMAVEn_dir, "integer_cuts_test")
    KS.MAPs_part_clashes(cut_file, 1.0)
else:
    cat1 = sys.argv[1]
    cat2 = sys.argv[2]
    cut_file = sys.argv[3]
    cutoff = float(sys.argv[4])
    KS.MAPs_part_category_clashes(cat1, cat2, cut_file, cutoff)
