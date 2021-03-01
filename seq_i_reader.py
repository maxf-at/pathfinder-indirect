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

import pathfinder_i
import pathfinder
from helper import print_moves

import pandas as pd


# filename = f'./sample_seqs/indirect_60.csv'
# filename = f'./sample_seqs/indirect_100.csv'
filename = f'./sample_seqs/indirect_ea_60.csv'
# filename = f'./sample_seqs/indirect_ea_100.csv'


df = pd.read_csv(filename)


Verbose = False
Verbose = True
search_width = 500
indirect_iterations = 2

Debug = False


for index, row in df.iterrows():
    
    add_moves = []

    # 60.7

    # 100.7:
    # add_moves = [(73,82), (-57,-94)]

    # 100.15 requires high search width


    if index!=10:
        continue
    # if index<14:
    #     continue


    sequence = row.sequence
    s1 = row.s1
    s2 = row.s2
    
    print (sequence)
    print (s1)
    print (s2)


    # launch pathfinder_i
    pathfinder_i.find_path(sequence, s1, s2, indirect_iterations=indirect_iterations, add_moves=add_moves, search_width=search_width, Debug=Debug, Verbose=Verbose)
    
    print (index, "indirect result:", row.max_en_i, "direct:", row.max_en_d)
    

    filename_dir = filename.split(".csv")[0]
    path_filename = f'{filename_dir}/{sequence}.txt'



    print ("\nexisting result:")
    # parse existing eapath / tabupath file
    if os.path.exists(path_filename):
        with open(path_filename, "r") as f:
            
            moves = []
            sequence = ""
            s1 = ""
            s2 = ""
            max_en = float("-inf")
            for line in f:
                line = line.strip()
                line = line.split(" ", 1)
                if sequence == "":
                    sequence = line[0]
                    continue
                if s1 == "":
                    s1 = line[0]
                if len(line)>1:
                    s2 = line[0]
                    rest = line[1]                    
                    rest = rest.replace('(', '').replace(')', '')
                    rest = rest.replace('[', '').replace(']', '')
                    rest = rest.replace('{', '').replace('}', '')
                    rest = rest.replace(',', '')
                    rest = rest.split()
                    if len(rest)!=3:
                        # print (line)
                        continue
                    i, j, en = rest
                    i = int(i)
                    j = int(j)
                    en = float(en)
                    if en > max_en:
                        max_en = round(en, 2)
                    moves.append((i, j, en))
  
            print_moves(sequence, s1, s2, moves)


        

    sequence = row.sequence
    s1 = row.s1
    s2 = row.s2

    print ("\nregular findpath result:")
    pathfinder_result = pathfinder.pathfinder(sequence, s1, s2, search_width=search_width, verbose=True)
    # print ("regular findpath result:", pathfinder_result.max_en)


    # if index>0:
    #     break
    

 

