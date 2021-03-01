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

import pathfinder
import helper
import merge_composition

import pandas as pd

# generate X random sequences ~ x nt for merge join

# target_count = 200
# min_bp_dist = 80  # 50  / 80
# max_bp_dist = 160 # 120 / 160
# x = 300 # length 300

# target_count = 100
# min_bp_dist = 200
# max_bp_dist = 300
# x = 500


# target_count = 1000
# min_bp_dist = 100
# max_bp_dist = 200
# x = 400


# target_count = 1000
# min_bp_dist = 80
# max_bp_dist = 160
# x = 300


# target_count = 1000
# min_bp_dist = 60
# max_bp_dist = 120
# x = 200

# target_count = 500
# min_bp_dist = 40
# max_bp_dist = 80
# x = 150

# target_count = 400
# min_bp_dist = 10
# max_bp_dist = 30
# x = 60

# indirect generation

# target_count = 100
# min_bp_dist = 10
# max_bp_dist = 20
# x = 70
target_count = 500
min_bp_dist = 30
max_bp_dist = 40
x = 80
# target_count = 100
# min_bp_dist = 20
# max_bp_dist = 30
# x = 100
# target_count = 100
# min_bp_dist = 40
# max_bp_dist = 60
# x = 150

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
    if bp_dist > max_bp_dist:
        continue
    l = merge_composition.merge_check(sequence, s1, s2, Debug=False)

    # ignore the case where sections provide no benefit
    if len(l) == 1 or l==[0,x]:
        continue


    # for the time being, only consider [x, [i,j], y]
    # if len(l[1]) > 2 or len(l)>3:
    #     continue

    # ignore [x, [i,j], y]
    # if len(l[1]) <= 2 and len(l)<=3:
    #     continue

    

    # if len(l[1])==2:
    #     continue


    print (len(s_list), len(l), bp_dist, l)

    search_width=100


    init_verbose = True

    init_results = pathfinder.pathfinder(sequence, s1, s2, search_width=search_width, verbose=init_verbose)



    s_list.append((sequence, s1, s2, l))





    
# run pathfinder, save results into csv

for sequence, s1, s2, l in s_list:

    # p,e = merge_dp.recursive_findpath(sequence, s1, s2, sections=l, search_width=5000, Verbose=Verbose)
    # print (e)

    sE_250 = 0
    sE_250 = pathfinder.pathfinder(sequence, s1, s2, search_width=250, verbose=False).max_en
    sE_250 = round(sE_250,2)
    print ("250 done", sE_250)

    sE_500 = pathfinder.pathfinder(sequence, s1, s2, search_width=500, verbose=False).max_en
    sE_500 = round(sE_500,2)
    # sE_500 = 0
    print ("500 done", sE_500)

    sE_5000 = pathfinder.pathfinder(sequence, s1, s2, search_width=5000, verbose=False).max_en
    sE_5000 = round(sE_5000,2)
    # sE_5000 = 0
    print ("5000 done", sE_5000)

    # sE_10000 = pathfinder.pathfinder(sequence, s1, s2, search_width=10000, verbose=False).sE
    # print ("10000 done", sE_10000)
    # e = round(e,2)
    
    # if e > sE_5000:
        # print ((sequence,s1,s2,sE_500, sE_5000, sE_10000))
        # s_list.append((sequence,s1,s2,sE_500, sE_5000, sE_10000))

    # print ((sequence,s1,s2,sE_250,sE_500, sE_5000))
    results.append((sequence,s1,s2,sE_250,sE_500, sE_5000))


df = pd.DataFrame(results, columns=['sequence','s1', 's2', 'sE_250', 'sE_500', "sE_5000"]).set_index('sequence') 
# filename = f'../sample_seqs/{len(df)}x_{x}-{min_bp_dist}.csv'
filename = f'./sample_seqs/indirect_input_80.csv'
df.to_csv(filename)
# df.to_csv(filename, mode='a', header=False)

# print (df)
