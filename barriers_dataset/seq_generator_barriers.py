#!/usr/bin/env python3
# coding: utf-8

import RNA
import helper
import pandas as pd
import pathfinder

# generate X random sequences ~ x nt for merge join


# indirect generation
# target_count = 100
# min_bp_dist = 10
# max_bp_dist = 20
# x = 70
# target_count = 100
# min_bp_dist = 20
# max_bp_dist = 30
# x = 100
# target_count = 100
# min_bp_dist = 40
# max_bp_dist = 60
# x = 150

target_count = 100
min_bp_dist = 10
max_bp_dist = float("inf")
en_range_cutoff = 10
x = 150 # seq length

# target_count = 100
# min_bp_dist = 10
# max_bp_dist = float("inf")
# en_range_cutoff = 10
# x = 100 # seq length




results = []

# check confirmations for merging

def detect_local_minimum(fc, structure):
    # perform gradient walk from sample to determine direct local minimum
    pt = RNA.IntVector(RNA.ptable(structure))
    fc.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
    return RNA.db_from_ptable(list(pt))



useless_sections = 0

while (len(results) < target_count):
    
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



    sw = 2
    init_verbose = False

    init_result = pathfinder.pathfinder(sequence, s1, s2, search_width=bp_dist*sw, verbose=False)

    # compute minimum free energy (MFE) and corresponding structure
    ss, mfe = RNA.fold(sequence)
    
    # calculate the energy range with respect to the mfe for RNAsubopt
    min_energy = min(init_result.e1_en,  init_result.e2_en, mfe)
    en_range = init_result.max_en - mfe

    if en_range > en_range_cutoff:
        continue

    # print (sequence)
    # print (s1)
    # print (s2)
    # print (bp_dist)
    # print (ss, mfe, "MFE")
    # print (init_result.max_en, en_range, mfe,  init_result.e1_en,  init_result.e2_en)


    elements = target_count
    index = len(results)

    percent_complete = 100-(elements-index)/elements*100
    bar_length = 20
    filled_length = int(percent_complete/100*bar_length)
    rest = bar_length - filled_length
    bar = "█" * filled_length + '_' * rest
    print(
        f'\rComputing |{bar}| {percent_complete:.1f}% complete {index} ', end="\r")

    results.append((sequence, s1, s2, bp_dist, init_result.max_en))



print ()

df = pd.DataFrame(results, columns=['sequence','s1', 's2', "bp_dist", "findpath_direct"]).set_index('sequence') 

filename = f'./sample_seqs/{len(df)}x_{x}-{min_bp_dist}.csv'
# filename = f'../sample_seqs/local_min_200x1000_multiple_sections_min10.csv'

df.to_csv(filename)
print (df)
