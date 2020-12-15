#!/usr/bin/env python3
# coding: utf-8

import time
import subprocess
import os
import numpy as np
import pandas as pd

import glob
import os.path
from os import path

# import pathfinder
# import merge_cpp
# import RNA
# import helper

# import networkx as nx
# import matplotlib.pyplot as plt

# from helper import path_class

def parse_dir(directory, output_filename):



    files = directory + "/*txt"
    
    files = glob.glob(files)

    collect_data = []

    for input_file in files:

        max_en_i = -999
        max_en_d = -999

        s1 = ""
        s2 = ""
        sequence = ""

        with open(input_file, 'r') as f:
            for line in f:
                line = line.strip()
                print (line)
                if line[0] == "A" or line[0] == "U" or line[0] == "C" or line[0] == "G":
                    sequence = line
                elif s1 == "" and sequence != "" and line[0]!="S" and line[0]!="r": 
                    s1 = line.split()[0]
                elif s1 != "" and sequence != "" and line[0]!="S" and line[0]!="r":
                    s2 = line.split()[0]

                if line[0] == "S" and max_en_i == -999:
                    line = line.split()[1]
                    max_en_i = float(line)
                elif line[0] == "S" and max_en_i != -999:
                    line = line.split()[1]
                    max_en_d = float(line) 


        if max_en_i < max_en_d:
            collect_data.append((sequence, s1, s2, max_en_i, max_en_d))
            print (input_file, max_en_i, max_en_d)
            print (sequence)
            print (s1)
            print (s2)
            print ("~~~")
        

  
    df = pd.DataFrame(collect_data, columns=['sequence','s1', 's2', 'max_en_i', 'max_en_d']).set_index('sequence') 
    
    if path.exists(output_filename):
        df.to_csv(output_filename, mode='a', header=False)
    else:
        df.to_csv(output_filename)








if __name__ == '__main__':
  
    # directory = r"/mnt/c/Users/Max/Google Drive/arbeit/TABU/Output/100"
    directory = r"/mnt/c/Users/Max/Google Drive/arbeit/RNAEAPath/Output_100"
    
    output_filename = r"/mnt/c/Users/Max/Google Drive/pathfinder/sample_seqs/indirect_ea_100.csv"

    parse_dir(directory,output_filename)

