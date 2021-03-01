#!/usr/bin/env python3
# coding: utf-8

import time
import subprocess
import os
import numpy as np
import pandas as pd


import pathfinder
# import merge_cpp
import RNA
import helper

import networkx as nx
import matplotlib.pyplot as plt

from helper import path_class
from helper import print_moves

from sample_db import RNA_db



def save_moves(sequence, s1, s2, dataset, indirect_moves, description):

    output_file = "./sample_seqs/indirect_db.csv"
    direct_search_width=1000

    print ("saving:")

    direct_result = pathfinder.pathfinder(sequence, s1, s2, search_width=direct_search_width)
    direct_moves = [(i[0], i[1]) for i in direct_result.paths[0]]

    print ("direct path:")
    en_direct = print_moves(sequence, s1, s2, direct_moves)
    print ("indirect path:")
    en_indirect = print_moves(sequence, s1, s2, indirect_moves)

    # en_direct = print_moves(sequence, s1, s2, direct_moves)
    # en_indirect = print_moves(sequence, s1, s2, indirect_moves)

    # columns = ["sequence","s1","s2", "en_direct", "en_indirect", "direct_moves","indirect_moves","description"]
    # df = pd.DataFrame(columns=columns)
    # df.loc[0] = [sequence,s1,s2,en_direct, en_indirect, direct_moves,indirect_moves,description]

    # only save if indirect en is lower than direct en
    if en_indirect < en_direct:           
        rna_db = RNA_db(dbname='./sample_db/rna_samples.sqlite')  
        rna_db.insert_path(sequence, s1, s2, dataset, indirect_moves, description)

        # if os.path.isfile(output_file):
        #     df.to_csv(output_file, mode='a', header=False)
        # else:
        #     df.to_csv(output_file, mode='w', header='column_names')



def RNAtabu(sequence, s1, s2, dataset):


    cmd = "./get_barrier"

    
    path = "/mnt/c/Users/Max/TABU/"

    input_file = "/mnt/c/Users/Max/TABU/input.faa"
    


    with open(input_file, 'w') as f:
        f.write('>Input\n')
        f.write(sequence)
        f.write('\n')
        f.write(s1)
        f.write('\n')
        f.write(s2)
        f.write('\n')

    
    # cmd = path + cmd
    cmd = cmd + " -f ./input.faa 100 10 70 16"
  
    # save current working path 
    old_current_path = os.path.abspath(os.getcwd())
    os.chdir(path)

    print (cmd)

    output = subprocess.check_output(cmd, shell=True, encoding="utf8")
    
    # parse generated result file

    # print (output)
    output = output.split("\n")



    fc = RNA.fold_compound(sequence)

    
    

    # p_result = pathfinder.pathfinder(sequence, s1, s2, search_width=1000)
    # p_en = p_result.max_en

    
    max_en = float("-inf")
    all_max_en = []

    current_path=[]
    all_paths=[]

    # TABU path outputs a bunch of paths as output, we search for the one with
    # the lowest saddle energy and save this one, discard the rest.

    for line in output:
        if line=="": continue
        # print (line)
        if line[0]==sequence[0]:          
            
            current_path=[]            
            max_en = float("-inf")     

        if line[0]=="b":            
            # print (max_en, "done")
            all_paths.append(current_path)
            all_max_en.append(max_en)
        if line[0] == "(" or line[0] == ")" or line[0] == ".":
            en = round(fc.eval_structure(line),2)
            current_path.append((line, en))
            # print (en)
            if en>max_en:
                max_en = en


    # only print/save the best path
    best_path = np.argmin(all_max_en)
    max_en = min(all_max_en)
    last_structure = s1
    

    moves = []

    for structure, en in all_paths[best_path]:
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
    save_moves(sequence, s1, s2, dataset, moves, "TABU")

    


    # # old code below to save output as txt file...
    #
    # output_file = f"./Output/{sequence}.txt"
    # e1 = fc.eval_structure(last_structure)
    # with open(output_file, 'w') as f:
    #     f. write(f'{sequence}\n')
    #     for structure, en in all_paths[best_path]:
    #         # calculate ptable move between structures
    #         diff = [i for i in range(len(structure)) if structure[i] != last_structure[i]]
    #         if last_structure != structure:
    #             if last_structure[diff[0]] != ".":
    #                 diff = (-diff[0]-1, -diff[1]-1)
    #             else:
    #                 diff = (diff[0]+1, diff[1]+1)
    #         else:
    #             diff = (0,0)
    #         moves.append(diff)
    #         last_structure = structure
    #         print (f'{structure} ({diff[0]:3}, {diff[1]:3}) {en:4.6}')
    #         f. write(f'{structure} ({diff[0]:3}, {diff[1]:3}) {en:4.6}\n')
    #     print (f'S: {max_en:6.6} kcal/mol, B:{max_en-e1:6.6} kcal/mol')
    #     print(f'regular findpath:')
    #     print(f'S: {p_en:6.6} kcal/mol, B:{p_en-e1:6.6} kcal/mol')
    #     f.write(f'S: {max_en:6.6} kcal/mol, B:{max_en-e1:6.6} kcal/mol\n')
    #     f.write(f'regular findpath:\n')
    #     f.write(f'S: {p_en:6.6} kcal/mol, B:{max_en-e1:6.6} kcal/mol\n')




    
    
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






    # print (result)










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

    # sequence = "CCCAAAGGGCCCAAAGGGCCCAAAGGGCCCAAAGGG"
    # s1 =       "(((...)))(((...)))(((...)))(((...)))"
    # s2 =       "(((............)))(((............)))"

    # sequence = "CGCGACGGCUACGCGACGGCAAUGCCGUUGCGAAGCCGUCGCGAUC"
    # s1 =       "(((((((((..............))))))))).............."
    # s2 =       "...........(((((((((..............)))))))))..."

        # rb2
    # sequence = "GGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCGUAAACUCUUGAUUAUGAAGUCUGUCGCUUUAUCCGAAAUUUUAUAAAGAGAAGACUCAUGAAU"
    # s1       = "............((((((.........))))))........((((((.......)))))).((((((.((((.((.((((((..........)))))).)))))))))))).." 
    # s2       = "(((((((((...((((((.........))))))........((((((.......))))))..))))))))).........................................." 

    # rb3
    # sequence = "CUCUUAUCAAGAGAGGUGGAGGGACUGGCCCGAUGAAACCCGGCAACCAGCCUUAGGGCAUGGUGCCAAUUCCUGCAGCGGUUUCGCGUUGAAAGAUGAGAGAUUCUUGUAGUCUCUUCUUUUAGCGAAGGGACUUUUUUU"
    # s1       = "((((((((....(.(((...(((.....)))......))))(((..((((((....))).))).)))........(((((......)))))...)))))))).........((((((((.......))))))))......."
    # s2       = "............(.(((...(((.....)))......))))(((..((((((....))).))).)))........(((((......)))))(((((.((((...........)))).)))))..................."

    # rb5 -26.40 kcal/mol vs -23.60 kcal/mol direct
    # sequence = "GGGAAUAUAAUAGGAACACUCAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGACUAUGGGUGAGCAAUGGAACCGCACGUGUACGGUUUUUUGUGAUAUCAGCAUUGCUUGCUCUUUAUUUGAGCGGGCAAUGCUUUUUUUAUUCUCAUAACGGAGGUAGACAGGAUGGAUCCACUGA"
    # s1       = "................((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))........(((((........)))))............((((((((((((((.......))))))))))))))..........................................."
    # s2       = ".....................(((...(.(((((.......))))).)........((((((.......))))))..)))((((((((((((.(((((........)))))..............))))))))))))................................................................"

    # dsrA
    # sequence = "ACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUAAGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUU"
    # s1       = "..(((((((.....)))))))...(((((((((((((.......))))))).)))))).((((((((((.....))))))))))."
    # s2       = "..(((((((.....)))))))................(...(((((....)))))...)((((((((((.....))))))))))."



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


    # RNAtabu(sequence, s1, s2)
    
 
    # 

    # # filename = r"../sample_seqs/indirect_input_100.csv"
    

    # dataset = "indirect_60.csv"
    # dataset = "indirect_100.csv"
    # dataset = "indirect_ea_60.csv"
    dataset = "indirect_input_80.csv"

    filename = f"./sample_seqs/{dataset}"

    df = pd.read_csv(filename)

    for index, row in df.iterrows():
        # if index>1:
        #     break
        sequence = row.sequence
        s1 = row.s1
        s2 = row.s2
        print (sequence)
        print (s1)
        print (s2)
        RNAtabu(sequence, s1, s2, dataset)
        # break



