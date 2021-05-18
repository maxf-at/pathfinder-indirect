#!/usr/bin/env python3
# coding: utf-8

import RNA
import numpy as np

import subprocess
import matplotlib.pyplot as plt
from PIL import Image
from IPython.display import SVG, display
from collections import Counter
from collections import defaultdict

import difflib
import sys
import os
import random
import string
import time

import findpath
import helper
import merge_composition

import pandas as pd

import xplorer_new

# indirect generation

# target_count = 100
# min_bp_dist = 10
# max_bp_dist = 20
# x = 70
# target_count = 500
# min_bp_dist = 30
# max_bp_dist = 40
# x = 80
# target_count = 100
# min_bp_dist = 20
# max_bp_dist = 30
# x = 100

target_count = 100
min_bp_dist = 30
x = 120

Verbose = False

s_list = []
results = []


def detect_local_minimum(fc, structure):
    # perform gradient walk from sample to determine direct local minimum
    pt = RNA.IntVector(RNA.ptable(structure))
    fc.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
    return RNA.db_from_ptable(list(pt))

# check confirmations for merging


while (len(s_list) < target_count):

    # RNAsubopt wrapper
    sequence, s1, s2 = helper.generate_structures(x)

    # overwrite s1 s2 to local minimum
    fc = RNA.fold_compound(sequence)
    s1 = detect_local_minimum(fc, s1)
    s2 = detect_local_minimum(fc, s2)

    bp_dist = RNA.bp_distance(s1, s2)

    if bp_dist < min_bp_dist:
        continue
    # if bp_dist > max_bp_dist:
    #     continue
    l = merge_composition.merge_check(sequence, s1, s2, Debug=False)

    print(f"sequence = '{sequence}'")
    print(f"s1       = '{s1}'")
    print(f"s2       = '{s2}'")

    

    # swm = 100
    # init_verbose = False
    # init_result = findpath.init_single_findpath(sequence, s1, s2, swm, True)
    # init_result = round(init_result/100.0,2)
    # print(len(s_list), len(l), bp_dist, l, "init result:", init_result)

    print(len(s_list), len(l), bp_dist, l)

    e1, e2, intermediate_structures = xplorer_new.init(sequence, s1, s2, 100)

    print(e1, e2, len(intermediate_structures))

    if e1 <= e2:
        print ("ignore")
        continue

    s_list.append((sequence, s1, s2, e1, e2, intermediate_structures))
    # break

# run pathfinder, save results into csv

# for sequence, s1, s2, l in s_list:
#     results.append((sequence, s1, s2))


df = pd.DataFrame(s_list, columns=['sequence', 's1', 's2', 'e1', 'e2', 'intermediate_structures']).set_index('sequence')

filename = f'./sample_seqs/indirect_new_120.csv'
# df.to_csv(filename)

print(df)
