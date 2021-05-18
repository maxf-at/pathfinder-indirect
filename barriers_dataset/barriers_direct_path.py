#!/usr/bin/env python3
# coding: utf-8

import sys
import random
import string
import RNA
import numpy as np
import pandas as pd
import os
import subprocess
from dataclasses import dataclass, field

import pathfinder




def barriers_indirect(index, sequence, s1, s2):

    bp_dist = RNA.bp_distance(s1,s2)
    Verbose = True
    # Verbose = False
    fp_result = pathfinder.pathfinder(sequence, s1, s2, search_width=bp_dist*20, verbose=Verbose)
    # barrier = result.max_en - result.e1_en

    result = fp_result

    pt1 = RNA.ptable_from_string(s1)
    pt2 = RNA.ptable_from_string(s2)

    move_list = []
    for i in range(1, pt1[0]):
        if (pt1[i] != pt2[i]):
            if (i < pt1[i]):
                move_list.append((-i, -pt1[i]))

            
            if (i < pt2[i]): 
                # insert this pair
                move_list.append((i, pt2[i]))


    # compute minimum free energy (MFE) and corresponding structure
    ss, mfe = RNA.fold(sequence)

    print (ss, mfe, "MFE")

    # calculate the max energy range with respect to the mfe
    min_energy = min(result.e1_en,  result.e2_en, mfe)
    barrier = result.max_en - mfe


    # barrier = 10.0

    print (result.max_en, barrier, mfe,  result.e1_en,  result.e2_en)

    index = f'{index}'

    if not os.path.exists(index):
        os.mkdir(index)
    os.chdir(index)


    if not os.path.exists('RNAsubopt_output_direct.txt') or os.stat("RNAsubopt_output_direct.txt").st_size == 0:
        # RNAsubopt output does not exist

        with open('barriers_input.txt', 'w') as f:
            f.write(f'{sequence}\n')

        cmd = f"RNAsubopt -d2 -s -e {barrier} < barriers_input.txt > RNAsubopt_output.txt"
        result = subprocess.check_output(cmd, shell=True, encoding="utf8")
        # print (result)

        print ("subopt done")

        cmd = f'../../../remove_indirect "{s1}" "{s2}"'
        result = subprocess.check_output(cmd, shell=True, encoding="utf8")

        print (result)

        # # delete indirect moves from subopt output
        # with open('RNAsubopt_output_direct.txt', 'w') as f1:
        #     with open('RNAsubopt_output.txt', 'r') as f2:
        #         for line_index, line in enumerate(f2):
                
        #             if line_index == 0:
        #                 f1.write(line)
        #                 continue

        #             l = line.strip()
        #             current_s = l.split()[0]
        #             # print(current_s, l)
        #             pt = RNA.ptable_from_string(current_s)

        #             indirect = False
        #             for i in range(1, pt[0]):
        #                 if (i < pt[i]): # ()                        
        #                     j = pt[i]
        #                     # if we have a basepair which is not present in s1 or s2 -> ignore
        #                     if (pt[i]!=pt1[i] or pt[j]!=pt1[j]) and (pt[i]!=pt2[i] or pt[j]!=pt2[j]):
        #                         indirect = True
                            
                            
        #                 if (pt[i] == 0 and pt1[i] != 0 and pt1[i]==pt2[i]):
        #                     indirect = True
        #                     # print (current_s, "i", i, pt[i], pt1[i], pt2[i])

        #             if (not indirect):
        #                 f1.write(line)

        #             # if line_index > 10:
        #             #     break

        # subopt done, delete old files
    if os.path.exists('RNAsubopt_output.txt'):
        os.remove('RNAsubopt_output.txt')

    all_files = os.listdir("./")

    for item in all_files:
        if item.startswith("path."):
            print ("delete", item)
            os.remove(item)


    # cmd = f"barriers --rates -G RNA --max 50 --minh 0.1 < RNAsubopt_output.txt > barriers_output.txt"
    # cmd = f"barriers --rates -G RNA --max 50 --minh 0.1 < RNAsubopt_output_direct.txt > barriers_output.txt"
    cmd = f"barriers --rates -G RNA --max 400 < RNAsubopt_output_direct.txt > barriers_output.txt"
    # cmd = f"barriers --rates -G RNA < RNAsubopt_output_direct.txt > barriers_output.txt"
    result = subprocess.check_output(cmd, shell=True, encoding="utf8")
    # print (result)

    id1 = False
    id1_s = ""
    id1_bp = 9999
    id2 = False
    id2_s = ""
    id2_bp = 9999

    with open('barriers_output.txt', 'r') as f:
        for i, line in enumerate(f):
            line = line.strip().split()
            if line[0][0] in ["S", "A", "C", "G", "T", "U"]:
                continue
            structure = line[1]
            bp1 = RNA.bp_distance(structure, s1)
            if bp1<id1_bp:
                id1 = int(line[0])
                id1_s = line[1]
                id1_bp = bp1     

            bp2 = RNA.bp_distance(structure, s2)  
            if bp2<id2_bp:
                id2 = int(line[0])
                id2_s = line[1]
                id2_bp = bp2

            # break

    print ("Structures: s1:", id1, id1_s, RNA.bp_distance(id1_s, s1))
    print ("Structures: s2:", id2, id2_s, RNA.bp_distance(id2_s, s2))
    print ("distance:", RNA.bp_distance(id1_s, id2_s))

    # cmd = f"barriers -P {id1}={id2} --rates -G RNA --max 50 --minh 0.1  < RNAsubopt_output.txt > barriers_output.txt"
    # cmd = f"barriers -P {id1}={id2} --rates -G RNA --max 50 --minh 0.1  < RNAsubopt_output_direct.txt > barriers_output.txt"
    cmd = f"barriers -P {id1}={id2} --rates -G RNA --max 400 < RNAsubopt_output_direct.txt > barriers_output.txt"
    # cmd = f"barriers -P {id1}={id2} --rates -G RNA < RNAsubopt_output_direct.txt > barriers_output.txt"


    print (cmd)
    result = subprocess.check_output(cmd, shell=True, encoding="utf8")
    print (result)

    os.chdir("../")

    return fp_result.max_en
    

def parse_barriers_output(index, sequence):

    index = f'{index}'

    moves = [(0,0)]

    structures = []

    s1 = ""
    s2 = ""
    last_s = ""
    last_ptable = False

    for subdir, dirs, files in os.walk(index):
        for file in files:
            if file[0:4] != "path":
                continue
            
            # open the path_xxx file
            with open(f"./{index}/{file}", "r") as f:                               

                for line in f:
                    line = line.strip().split()
                    s = line[0]
                    structures.append(s)

                    if not last_ptable:
                        s1 = s
                        last_s = s
                        last_ptable = RNA.ptable_from_string(s)
                    
                    current_ptable = RNA.ptable_from_string(s)
                    s2 = s
                         
                    bp_diff = RNA.bp_distance(last_s, s)
                    print (s, moves[-1], bp_diff)

                    # sometimes, barriers creates shift moves... lets fix it.

                    if bp_diff > 1:
                        fix_moves = pathfinder.pathfinder(sequence, last_s, s, search_width=bp_diff*20, verbose=False).paths
                        fix_moves = fix_moves[0][1:] # delete (0,0) move
                        fix_moves = [(x[0], x[1]) for x in fix_moves]
                        
                        for x,y in fix_moves:
                            moves.append((x, y))                        
                        print (fix_moves)

                    else:
                        # use moves we found with barriers
                        for i in range(1,len(sequence)):
                            if last_ptable[i] != current_ptable[i]:                            
                                if current_ptable[i] == 0 and i<last_ptable[i]: # del bp
                                    moves.append((-i, -last_ptable[i]))                            
                                if last_ptable[i] == 0 and i<current_ptable[i]: # add
                                    moves.append((i, current_ptable[i]))

                    last_ptable = current_ptable
                    last_s = s

    print ("found structures in path:", (len(structures), len(moves)))

    return s1, s2, structures, moves


def print_moves(sequence, s1, s2, moves, move_color='\033[93m', structures=False, Verbose = True, exclude=None, include=None):

    """
    print a folding path with colour coding
    
    moves have to contain i,j or optionally i,j,en
    e.g. [(0, 0, -15.4), (-2, -24, -13.2), (-3, -23, -11.6)]

    without verbosity, this just returns max_en
    """

    # print (moves)

    # from stackexchange...
    class c:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        CYAN = '\033[96m'
        GREEN = '\033[92m'
        YELLOW = '\033[93m'
        RED = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'

    if Verbose: print(sequence)
    s = s1

    move_color='\033[92m'

    fc = RNA.fold_compound(sequence)
    e1 = en = round(fc.eval_structure(s), 2)
    max_en = float("-inf")

    output_rows = []
    moves_i_j = [(x[0],x[1]) for x in moves]

    indirect_moves = []
    direct_moves = []

    # preprocessing - generate strings & energies if required
    # for a in range(len(structures)):
    for a in range(len(moves)):

        move = moves[a]
        en = False
        if len(move) == 2:
            i, j = move
        if len(move) == 3:
            i, j, en = move
        if structures:
            s = structures[a]
        else:
            if i > 0:
                s = s[:i-1] + "(" + s[i:j-1] + ")" + s[j:]
            if i < 0:
                s = s[:-i-1] + "." + s[-i:-j-1] + "." + s[-j:]
        # print("current s:", s)
        if not en:
            en = round(fc.eval_structure(s), 2)
        e2 = en
        if en > max_en:
            max_en = en
        output_rows.append((s, i, j, en))


    for s, i, j, en in output_rows:

        # print initial row with move (0,0)
        if i == 0:
            info = f'{move_color}[{i:4}, {j:4} ]{c.ENDC} {en:6.2f}'
            if Verbose: print(f"{s} {info}")
            continue

        pos_i = abs(i)-1
        pos_j = abs(j)-1

        # if a move has an identical inverse copy (1,2) <-> (-1,-2)
        # it is automatically an indirect move - these are colored in red 
        if (-i, -j) in moves_i_j: # indirect move
            if i>0:
                indirect_moves.append((i,j))
            colored_s = s[0:pos_i] + c.RED + c.BOLD + s[pos_i] + c.ENDC +\
                s[pos_i+1:pos_j] + c.RED + c.BOLD + \
                s[pos_j] + c.ENDC + s[pos_j+1:]
            info = f'{c.RED}[{i:4}, {j:4} ]{c.ENDC}'
        else:  # direct move
            colored_s = s[0:pos_i] + move_color + c.BOLD + s[pos_i] + c.ENDC +\
                s[pos_i+1:pos_j] + move_color + c.BOLD + \
                s[pos_j] + c.ENDC + s[pos_j+1:]
            info = f'{move_color}[{i:4}, {j:4} ]{c.ENDC}'

        if en == max_en:
            info += f' {c.RED}{c.BOLD}{en:6.2f}{c.ENDC}'
        else:
            info += f' {en:6.2f}'

        if Verbose: print(f"{colored_s} {info}")

    barrier = max_en - e1
    if Verbose: print(
        f"S: {max_en:6.2f} kcal/mol | B: {barrier:6.2f} kcal/mol | E[start]:{e1:6.2f} E[end]:{e2:6.2f}")
    
    return indirect_moves, output_rows, max_en





if __name__ == "__main__":

    # Tabu paper toy example
    # sequence = "CGCGACGGCUACGCGACGGCAAUGCCGUUGCGAAGCCGUCGCGAUC"
    # s1 = "(((((((((..............))))))))).............."
    # s2 = "...........(((((((((..............)))))))))..."

    # rna2dfold example
    # sequence = "GGGCGCGGUUCGCCCUCCGCUAAAUGCGGAAGAUAAAUUGUGUCU"
    # s1 = "(((((.....)))))(((((.....)))))(((((.....)))))"
    # s2 = "((((((((((.....(((((.....))))).....))))))))))"

    # sequence = "AAGAAGACCUCAAUCGAAUCACGGGCAAGUCCGACGAGGAACGCCUAGGCGAGGUGAUCGGCCCGAUCUUAAUGUAGGAU"
    # s1 = ".......(((...........(((((...(((.....))).(((((.....)))))....)))))..........))).."
    # s2 = "..((.((......))...)).(((((..((((..(((((....)))...))..).)))..)))))((((((...))))))"

    sequence = "CUUCCAGAGUUUCGUUGAACCGGCAAACGGCCACUAAGGAUAGCGCCUCCCGCUUGUCCCUGGACCGGCGUGCAUCGUAU"
    s1 = "..(((((.((((....)))).(((.....))).....((((((((.....)).)))))))))))...(((.....))).."
    s2 = "......((...)).........(((..((((((....((((((((.....)))).)))).))).)))...)))......."

    # sequence = "CUGGGUCGUGGUGGCCUUUUAGAUACGAUUCACGAACGUAGCACGUUUCGGUCUCCGGAGACGCAAUGAUCUCGAGGGUA"
    # s1 = ".(((((((((.(((....)))..))))))))).((.(((.((..(((((((...))))).)))).))).))........."
    # s2 = ".(((((((((.............)))))))))(((((((...))))).)).((((.(((..........))).))))..."

    # sequence = "GGAAGCCGGCGAGGCAGUACCAUUAUAUAGUUUGUCUUCCAAGAAUGGGUACGACCGCGGGACCGUUCGGUUAUCGUCUG"
    # s1 =       "((.((((((((.((..((((((((....((.....))......))))).)))..)).)).......)))))).))....."
    # s2 =       "((.((((((.(((((((.((.........)))))))))))..((((((...((....))...)))))))))).))....."

    # sequence = 'UAAAAUGAUCACGGUUUCAGCUUUGGACAGGGCGUUCCACUAAACUCCUGGUGACAUAGAUAUAUUGGAUUGCAACUACUCGUCGGUCCGGUUGGCGUUC'
    # s1       = '....(((.((((...............(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........'
    # s2       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'

    # sequence = "ACUUAACGUCUCCCUUAAUGAAGUGUUCCUCUAUAGUAACUAAGGAAAUACCGUGCGUCAUACCUGGAGUUCGGGGUACU"
    # s1 = "......((..((.......))..))(((((............)))))((.(((.((.(((....))).)).))).))..."
    # s2 = "((((.....((((....((((.((((((((...(((...))))))))......))).))))....))))....))))..."

    # sequence = 'UUGGGAAUCCAUCUUCGAGACAGUGUACAUAGAUACCAAUGAUCUAGAAACAAUUAUCAGCACUCCCGCUUCGUCGGGGAACCGGACGUAGUCUUAGAGAACAAUUCUUUCUACAAUUCGAAAAGCUACUACCGUGGUGGAAUCUAGAGA'
    # s1       = '...(((.((((((..((......(((...(((((.......)))))...)))...........(((((......)))))........(((((..((((((.......))))))...........)))))...)).)))))).))).....'
    # s2       = '..((((((...((((.(((((..(((...(((((.......)))))...)))......(((......))).((((((....)))).))..))))).))))...))))))((((..((((.....(((((....))))).)))).))))..'

    # sequence = 'AACAUUUAUAGGCGUAUGACGGGGUUUCACAGCUAGAUCCAUGAUGUAAUUGUGAUAGUCUACGUAUCAACCGUACCCCCGGGGAAAUCCUGAAAUGUAGCGCGUCACGACUAAAGUAAUUCAGUAUCGUCGGUGGUAGGCGUCUCAAAU'
    # s1       = '.((((((.((((.......(((((..((((((.((.(((...))).)).)))))).....((((.......)))).))))).......))))))))))((.(((((((.(((...(((......)).)...))).)).))))))).....'
    # s2       = '.((((((.((((.......(((((..((((((.((.(((...))).)).)))))).....((((.......)))).))))).......)))))))))).((((.(((((((....(((......))).))).))))...)))).......'

    # sequence = 'GAGCGCAGUGUACGCGUAUCAACCAAAGGCGGCCCAGUGGUUGCGCCGAGUACGCCAACGACAGACAGUCCCCGUUCGUCCGAUUAUCCAAGAACUAAUAUUAGAUAUGAGGUCGUAACAACCUCGCGUGAUAUCUCUAUUCAUGGUCAG'
    # s1       = '(((((..((((((.((...((((((..((....))..))))))...)).))))))....(((.....)))..)))))....((((((....(((...((((((....((((((.......))))))..)))))).....)))))))))..'
    # s2       = '(((((..((((((.((...((((((..((....))..))))))...)).))))))....(((.....)))..)))))................(((((((..(((((((((((.......))))).....)))))).)))...))))...'

    # #21
    # sequence = 'UCCGACAUUAAGACACACCAGGGUCUCGUAUCCCUAGGGUAAGGUACGCGCGGACCGGCCAAUCGGGUAUUGCUGCAAACUAUGGCAAUAGUGCAUAGGUUCAGACGAAGUACGGGUGGAUAUUUGUAGCCAGUAUGCUGGGUCUCCGGG'
    # s1       = '((((..........((.((((((........)))).))))..........))))((((...(((.(((((.((((((((((((.(((....)))))))(((((..((.....))..))))).))))))))....))))).)))..)))).'
    # s2       = '((((..........((.((((((........)))).))))..........))))((((....((((((((((((((((((((((((....)).)))))(((((..((.....))..))))).)))))))).))))).))))....)))).'

    # sequence = 'CUUCGGAGACCUAAGGAAUACGUCUUUCGAGCGGAAGGAUUCCUCGUUCAAUGCCAACUUAUUCUUGUUGGUUAAGGCAACGAAUCCCCUACGUCUCAACACGAAGUGGAUGACCAGGGAUGAAACCACGCUGGGCCAGCUUCCUCCCUC'
    # s1       = '....((((.....(((((....((((((.....)))))))))))........((((((........))))))...(((.....((((.((.(((......))).)).))))..((((.(.((....))).)))))))......))))...'
    # s2       = '((((((((((....(((....(((((((.....)))))))...((((((...((((((........))))))....).))))).))).....)))))....))))).........(((((.(((.....((((...))))))).))))).'

    # sequence = 'AGCGGAGGGGUCUAGUGGUCACAUCCUACGAGUUCGGGGCACCGAACAAAUCGAAGUCUACCGAGAUAGGGCGUUCGUUAGGCGGACGCUAGCCCAUCAUUUCUCAUGGUUCAUUGAUAC'
    # s1       = '....((.(((.((((((..(....(((((((((((((....)))))).............((.......))...)))).))).)..))))))))).))...((....))...........'
    # s2       = '...((.(((((((..((((.((.((......((((((....)))))).....)).))..))))))))..(((((((((...)))))))))..))).))...((..(((...))).))...'

    # barriers_indirect(sequence, s1, s2)#
    # s1, s2, structures, moves = parse_barriers_output(sequence)
    # print (moves)
    # print (f"sequence = '{sequence}'")
    # print (f"s1       = '{s1}'")
    # print (f"s2       = '{s2}'")
    # indirect_moves, max_en = print_moves(sequence, s1, s2, moves, structures=structures)
    # print (indirect_moves)


    
    # filename = f'./sample_seqs/100x_100-10.csv'
    # results_dir = "100x_100-10"
    # results_file = "100x_100-10.txt"

    filename = f'./sample_seqs/100x_150-10.csv'
    results_dir = "100x_150-10"
    results_file = "100x_150-10.txt"

    df = pd.read_csv(filename)

    new_dir = filename
    os.chdir("./results/")
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    os.chdir(results_dir)

    
    for index, row in df.iterrows():

        # if index!=58:
        #     continue
        # if index>10:
        #     break

        sequence = row.sequence
        s1 = row.s1
        s2 = row.s2

        print ()
        print ("index:", index)
        print (f"sequence = '{sequence}'")
        print (f"s1       = '{s1}'")
        print (f"s2       = '{s2}'")

        try:
            fp_max_en = barriers_indirect(index, sequence, s1, s2)
        except:
            print ("ignore", index)
            os.chdir("../")
            continue

        # s1 and s2 might be slightly altered...
        
        s1, s2, structures, moves = parse_barriers_output(index, sequence)
        
        try:
            # indirect_moves, move_tuples, barriers_max_en = print_moves(sequence, s1, s2, moves, structures=structures)
            indirect_moves, move_tuples, barriers_max_en = print_moves(sequence, s1, s2, moves)

        except:
            print ("index:", index, "invalid path")
            continue


        # if len(indirect_moves) > 0:
        #     print ("ignore, indirect moves found")
        #     continue


        bp_dist = RNA.bp_distance(s1, s2)

        print ("fp:", fp_max_en, "barriers:", barriers_max_en)

        results = [(index, sequence, s1, s2, bp_dist, fp_max_en, barriers_max_en, moves, indirect_moves)]
        df = pd.DataFrame(results, columns=['index', 'sequence', 's1', 's2', 'bp_dist', 'fp_max_en', 'barriers_max_en', 'moves', 'indirect_moves']).set_index('index') 

        # create new csv if required, otherwise append
        with open(results_file, 'a') as f:
            df.to_csv(f, header=f.tell()==0)  

