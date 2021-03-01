#!/usr/bin/env python3
# coding: utf-8

import time
import subprocess
import os
import numpy as np
import pandas as pd
import os.path


import pathfinder
import RNA
import helper

import networkx as nx
import matplotlib.pyplot as plt

from helper import path_class
from helper import print_moves
from sample_db import RNA_db


def save_moves(sequence, s1, s2, dataset, indirect_moves, description):

    direct_search_width=1000

    direct_result = pathfinder.pathfinder(sequence, s1, s2, search_width=direct_search_width)
    direct_moves = [(i[0], i[1]) for i in direct_result.paths[0]]

    print ("direct path:")
    en_direct = print_moves(sequence, s1, s2, direct_moves)
    print ("indirect path:")
    en_indirect = print_moves(sequence, s1, s2, indirect_moves)

    # only save if indirect en is lower than direct en
    if en_indirect < en_direct:           
        rna_db = RNA_db(dbname='./sample_db/rna_samples.sqlite')  
        rna_db.insert_path(sequence, s1, s2, dataset, indirect_moves, description)


def RNAeapath(sequence, s1, s2, dataset):


    cmd = "./RNAeapath"

    path = "/mnt/c/Users/Max/RNAEAPath/"

    input_file = "/mnt/c/Users/Max/RNAEAPath/input.faa"
    output_file = "./Output/log.txt"


    with open(input_file, 'w') as f:
        f.write('>Input\n')
        f.write(sequence)
        f.write('\n')
        f.write(s1)
        f.write('\n')
        f.write(s2)
        f.write('\n')

    
    # cmd = path + cmd
    # cmd = cmd + " -i ./input.faa  -o input.txt -d ./Output/ -c ./config.txt"
    cmd = cmd + " -i ./input.faa  -o input.txt -d ./Output/"
  
      # save current working path 
    old_current_path = os.path.abspath(os.getcwd())
    os.chdir(path)



    output_file = f"./Output/{sequence}.txt"
    if os.path.isfile(output_file):
        # print ("exists")
        return
        
    print (cmd)
    # return

    subprocess.check_output(cmd, shell=True, encoding="utf8")
    
    # parse generated result file
    output_file = "./Output/log.txt"

    result = []

    with open(output_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line!="" and line[0] == "G": # remove previous results
                result = []
            if line!="":
                result.append(line)


    fc = RNA.fold_compound(sequence)

    output_file = f"./Output/{sequence}.txt"
    last_structure = s1

    # p_result = pathfinder.pathfinder(sequence, s1, s2, search_width=1000)
    # p_en = p_result.max_en


    moves = []

    max_en = float("-inf")
    e1 = fc.eval_structure(last_structure)

    for line in result:        
        if line[0] == "(" or line[0] == ")" or line[0] == ".":
            line = line.split(" ")
            structure = line[0]
            en = fc.eval_structure(structure)
            if en>max_en:
                max_en = en
            # calculate ptable move between structures
            diff = [i for i in range(len(structure)) if structure[i] != last_structure[i]]
            if last_structure != structure:
                if last_structure[diff[0]] != ".":
                    diff = (-diff[0]-1, -diff[1]-1)
                else:
                    diff = (diff[0]+1, diff[1]+1)
            else:
                diff = (0,0)
            moves.append(diff)
            last_structure = structure

    os.chdir(old_current_path)
    save_moves(sequence, s1, s2, dataset, moves, "RNAeapath")

    
    # with open(output_file, 'w') as f:
    #     f. write(f'{sequence}\n')
        
    #     max_en = float("-inf")
    #     e1 = fc.eval_structure(last_structure)

    #     for line in result:        
    #         if line[0] == "(" or line[0] == ")" or line[0] == ".":
    #             line = line.split(" ")
    #             structure = line[0]
    #             en = fc.eval_structure(structure)
    #             if en>max_en:
    #                 max_en = en
    #             # calculate ptable move between structures
    #             diff = [i for i in range(len(structure)) if structure[i] != last_structure[i]]
    #             if last_structure != structure:
    #                 if last_structure[diff[0]] != ".":
    #                     diff = (-diff[0]-1, -diff[1]-1)
    #                 else:
    #                     diff = (diff[0]+1, diff[1]+1)
    #             else:
    #                 diff = (0,0)
    #             f. write(f'{structure} ({diff[0]:3}, {diff[1]:3}) {en:6.6}\n')
    #             print (structure, diff, en)
    #             last_structure = structure

    #     f. write(f'S: {max_en:6.6} kcal/mol, B:{max_en-e1:6.6} kcal/mol\n')
    #     f. write(f'regular findpath:\n')
    #     f. write(f'S: {p_en:6.6} kcal/mol, B:{max_en-e1:6.6} kcal/mol\n')
















if __name__ == '__main__':
    section = None

    
    # sbox_leader
    sequence = "AAAUGAUUCAACUUUUCAAAAAAUUAAUAACAUUUUCUCUUAUCGAGAGUUGGGCGAGGGAUUGGCCUUUUGACCCCAACAGCAACCGACCGUAAUACCAUUGUGAAAUGGGGCGCACUGCUUUUCGCGCCGAGACUGAUGUCUCAUAAGGCACGGUGCUAAUUCCAUCAGAUUGUGUCUGAGAGAUGAGAGAGGCAGUGUUUUACGUAGAAAAGCCUCUUUCUCUCAUGGGAAAGAGGCUUUUUGU"
    s1 =       "............................((((((((((((((((....((((((((((((......))))))..)))))).(((.(((.((......(((((....)))))(((((..........)))))(((((....)))))....))..))))))........((((((...))))))..))))))))))..)))))).......((((((((((((((((....)))))))))))))))).."
    s2 =       "............................(((((((((((((((((((.((((((.(((((.....)))))....))))))(((.((((.((......(((((....)))))(((((..........)))))(((((....)))))....))..)))))))..)))..((((((...))))))..)))))))))))..))))).......((((((((((((((((....)))))))))))))))).." 

    # s15 uninteresting
    sequence = "CUGGGAUCGCUGAAUUAGAGAUCGGCGUCCUUUCAUUCUAUAUACUUUGGAGUUUUAAAAUGUCUCUAAGUACU"
    s1 =       "..((((.((((((........))))))))))...........(((((.((((...........))))))))).."
    s2 =       ".((..(.((((((........))))))...)..)).......(((.((((((...........))))))))).."

    # ms2
    sequence = "GGGUGGGACCCCUUUCGGGGUCCUGCUCAACUUCCUGUCGAGCUAAUGCCAUUUUUAAUGUCUUUAGCGAGAC"
    s1 =       "(((..(((((((....)))))))..)))........(((..(((((.(.(((.....))).).)))))..)))"
    s2 =       "((((((((((((....))))))))))))........(((..(((((.((.(((...))))).)).)))..)))"

    sequence = "CCCAAAGGGCCCAAAGGGCCCAAAGGGCCCAAAGGG"
    s1 =       "(((...)))(((...)))(((...)))(((...)))"
    s2 =       "(((............)))(((............)))"

    sequence = "CGCGACGGCUACGCGACGGCAAUGCCGUUGCGAAGCCGUCGCGAUC"
    s1 =       "(((((((((..............))))))))).............."
    s2 =       "...........(((((((((..............)))))))))..."

    # >s15
    # CUGGGAUCGCUGAAUUAGAGAUCGGCGUCCUUUCAUUCUAUAUACUUUGGAGUUUUAAAAUGUCUCUAAGUACU
    # ..((((.((((((........))))))))))...........(((((.((((...........))))))))).. (-17.9)
    # ..((((.((((((........))))))))))...........((((..((((...........)))).)))).. (-14.3)
    # ..((((.((((((........))))))))))...........((((.(((((...........))))))))).. (-17.9)
    # ..((((.((((((........))))))))))...........(((..(((((...........))))).))).. (-13.9)
    # ..((((.((((((........))))))))))...........(((.((((((...........))))))))).. (-17.9)
    # ...(((.((((((........)))))))))............(((.((((((...........))))))))).. (-17.2)
    # ...((..((((((........)))))).))............(((.((((((...........))))))))).. (-12.4)
    # ...(...((((((........))))))..)............(((.((((((...........))))))))).. (-11.5)
    # .......((((((........))))))...............(((.((((((...........))))))))).. (-12.6)
    # ..(....((((((........))))))......)........(((.((((((...........))))))))).. (-10.8)
    # ..(..(.((((((........))))))...)..)........(((.((((((...........))))))))).. (-14.4)
    # .((..(.((((((........))))))...)..)).......(((.((((((...........))))))))).. (-15)


    # dataset = "test"
    # RNAeapath(sequence, s1, s2, dataset)

 

    # if parsing files...

    # 

    # filename = r"../sample_seqs/indirect_input.csv"
    # filename = r"../sample_seqs/indirect_input_100.csv"
    # filename = r"../sample_seqs/indirect_input_150.csv"

    # dataset = "indirect_60.csv"
    # dataset = "indirect_ea_60.csv"
    # dataset = "indirect_ea_100.csv"
    dataset = "indirect_input_80.csv"

    filename = f"./sample_seqs/{dataset}"

    df = pd.read_csv(filename)

    for index, row in df.iterrows():
        if index<338:
            continue

        sequence = row.sequence
        s1 = row.s1
        s2 = row.s2
        print (sequence)
        print (s1)#
        print (s2)
# 
        # if sequence == "UGUACGAUGUAUAUAGCUGGGAUAUGCAAUCACUUAUUCCCUGAGGAUGACUUCAUAUCUUUUACAGCGGCAGAAUUUGG":
            # print (index)



        RNAeapath(sequence, s1, s2, dataset)
        # break



