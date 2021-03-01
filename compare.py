#!/usr/bin/env python3
# coding: utf-8

# generic Python libraries
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from PIL import Image
from IPython.display import SVG, display
from collections import Counter
from collections import defaultdict
from dataclasses import dataclass
import itertools
import difflib
import sys
import os
import random
import string
import time
import logging

# coloredlogs
import coloredlogs

# custom
import RNA
from helper import p_to_s, print_moves
from sample_db import RNA_db
import pathfinder
import pathfinder_i



def find_sw_and_path(exceed=float("inf")):

    """
    compare search width indirect / direct findpath
    until optimal solution is found
    """

    max_en = float("inf")
    old_max_en = max_en

    search_width = 10
    Verbose = False
    Debug = False
    indirect_iterations = 1
    exceed_sw = 0

    while True:          
        # indirect_iterations = 2                

        max_en = pathfinder_i.find_path(sequence, s1, s2, indirect_iterations=indirect_iterations, add_moves=indirect_moves,
                                    search_width=search_width, Debug=Debug, Verbose=Verbose)
        # print ("iteration", search_width, max_en)
    
        if max_en < old_max_en:
            best_search_width = search_width

        if max_en < exceed:
            exceed_sw = search_width
            exceed=float("-inf") # trigger only once

        if search_width >= 500:
            # print (max_en)
            break

        search_width *= 1.414
        search_width = int(search_width)
        old_max_en = max_en

    return best_search_width, max_en, exceed_sw
    # print ("best sw:", )


def indirect_moves(sequence, s1, s2, moves):

    
    search_width = 500
    Verbose = False
    # Verbose = True
    Debug = False
    indirect_iterations = 1    

    # filter indirect moves from moves, use them as consecutive input
    indirect_moves = []    
    for i, j in moves:
        if (-i, -j) in moves and (-i, -j) not in indirect_moves and (i, j) != (0, 0):
            indirect_moves.append((i, j))
    # print ("input:", indirect_moves)

    max_en, best_path, used_indirect_moves = pathfinder_i.find_path(sequence, s1, s2, indirect_iterations=indirect_iterations, add_moves=indirect_moves,
                            search_width=search_width, Debug=Debug, Verbose=Verbose)

  
    # best_path_indirect_moves = []
    # # check which indirect moves were actually used
    # filtered = [(i[0], i[1]) for i in best_path]
    # for i, j in filtered:
    #     if (-i, -j) in filtered and (-i, -j) not in best_path_indirect_moves and (i, j) != (0, 0):
    #         best_path_indirect_moves.append((i, j))
    # best_path_indirect_moves.sort()

    print ("best case used indirect moves:", used_indirect_moves)
    # print_moves(sequence, s1, s2, best_path)



    # start at 10 and see from which point on we're using these indirect paths
    search_width = 10

    while search_width < 500:          
        # indirect_iterations = 2                            
            
        max_en, best_path, used_indirect_moves = pathfinder_i.find_path(sequence, s1, s2, indirect_iterations=indirect_iterations, add_moves=indirect_moves,
                                    search_width=search_width, Debug=Debug, Verbose=Verbose)



        # filtered = [(i[0], i[1]) for i in best_path]
        # current_indirect_moves = []
        # for i, j in filtered:
        #     if (i, j) in best_path_indirect_moves:
        #         current_indirect_moves.append((i, j))
        # current_indirect_moves.sort()

        print ("sw:", search_width, "indirect moves used:", used_indirect_moves, "en:", max_en)
        # print_moves(sequence, s1, s2, best_path)

        search_width *= 1.414
        search_width = int(search_width)







# 13 direct fp: 10 -5.7
# indirect fp: 10 -5.7 better after sw: 0 moves: [(-9, -23), (-8, -24), (-49, -74), (-7, -25), (6, 26), (5, 27), (-43, -80), (-7, -25), (-8, -24), (-9, 
# s1 s2 switch
#
# 18 direct fp: 548 -2.0
# indirect fp: 275 -0.9 better after sw: 0 moves: [(49, 60), (50, 59), (51, 58), (52, 57), (-3, -77), (44, 68), (45, 67), (46, 66), (47, 65), (-3, -77)]
# needs [(51, 58), (52, 57), (50, 59), (49, 60), (10, 17), (9, 18)]
#
# 27 direct fp: 388 -4.6
# indirect fp: 98 -4.1 better after sw: 0 moves: [(21, 39), (20, 40), (22, 38)]

# 33 direct fp: 548 1.1
# indirect fp: 388 1.3 better after sw: 0 moves: [(60, 71), (61, 70), (62, 69), (63, 68)]

# 48 direct fp: 70 -1.0
# indirect fp: 70 -1.0 better after sw: 0 moves: [(20, 53), (21, 52), (19, 54), (33, 41), (34, 40), (35, 39), (23, 52), (22, 53), (-26, -48)]

# 71 direct fp: 36 -4.6
# indirect fp: 36 -4.6 better after sw: 0 moves: [(-24, -48), (22, 48), (-25, -47), (23, 47), (-21, -49), (-20, -50), (-19, -51), (-18, -52), (-17, -53), (-17, -53), (-18, -52), (-19, -51), (-20, -50), (-21, -49)]


# 75 direct fp: 548 -5.8
# indirect fp: 138 -5.4 better after sw: 0 moves: [(28, 80), (29, 79), (30, 78), (31, 77), (-1, -25)]
# needs 
# add_moves = [(10, 18), (36, 75), (32, 79), (33, 78), (31, 80)]
# s1, s2 = s2, s1


# check when indirect paths are actually used?

if __name__ == '__main__':

    rna_db = RNA_db(dbname='./sample_db/rna_samples.sqlite')

    dataset = "indirect_input_80.csv"
    description = "RNAeapath"

    Verbose = True
    # Verbose = False


    for iteration, sequence, s1, s2, id, moves in rna_db.get_paths(dataset, description):

        if iteration <= 2:
            continue
        
        print ()
        print ("iteration", iteration)

        indirect_moves(sequence, s1, s2, moves)

        print (sequence)
        print (s1)
        print (s2)

        s1, s2 = s2, s1
        indirect_moves(sequence, s1, s2, moves)
        s1, s2 = s2, s1


        break





        # if iteration <= 10:
        #     continue

        # indirect_moves = []
        # # find saddle without any indirect moves

        # best_search_width, max_en, _ = find_sw_and_path()
        # print (iteration, "direct fp fwd:", best_search_width, max_en)

        # s1, s2 = s2, s1
        # best_search_width, max_en, _ = find_sw_and_path()
        # s1, s2 = s2, s1
        # print (iteration, "direct fp bwd:", best_search_width, max_en)

        # for i, j in moves:
        #     if (-i, -j) in moves and (-i, -j) not in indirect_moves and (i, j) != (0, 0):
        #         indirect_moves.append((i, j))

        # # indirect_moves = [(6,26),(5,27)]

        # # print(indirect_moves)
        # # db_saddle = print_moves(sequence, s1, s2, moves, Verbose=True)
        # db_saddle = print_moves(sequence, s1, s2, moves, Verbose=False)

        # # launch pathfinder_i with indirect moves from db data
        # # print(db_saddle)
        
        # best_search_width_i, max_en_i, exceed_sw = find_sw_and_path(exceed=max_en)
        # print (iteration, "indirect fp fwd:", best_search_width_i, max_en_i, "better after sw:", exceed_sw, "moves:", indirect_moves)

        # s1, s2 = s2, s1
        # best_search_width_i, max_en_i, exceed_sw = find_sw_and_path(exceed=max_en)
        # s1, s2 = s2, s1        
        # print (iteration, "indirect fp bwd:", best_search_width_i, max_en_i, "better after sw:", exceed_sw, "moves:", indirect_moves)


        # break