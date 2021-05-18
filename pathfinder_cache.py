#!/usr/bin/env python3
# coding: utf-8

# import sys
# import random
# import string
import RNA
import time
from dataclasses import dataclass

import matplotlib.pyplot as plt
from helper import path_class, print_moves

# @dataclass
# class return_values:
#     sE: int
#     en_min: float
#     s_pos: int
#     barrier: float
#     e1: float
#     e2: float
#     moves_str: list
#     moves_en: list
#     moves_pos: list
#     runtime: float


def pathfinder(sequence, s1, s2, verbose=False, output=False, search_width = None, section = None):

    # prune off sections which will currently not be regarded / later merged
    if section:
        if len(section)==2:
            start = section[0]
            end = section[1]
            s1 = '.'*start + s1[start:end] + '.'*(len(s1)-end)
            s2 = '.'*start + s2[start:end] + '.'*(len(s1)-end)
        if len(section)==4:
            # print ('section 4', section)
            start = section[0]
            mid1 = section[1]
            mid2 = section[2]
            end = section[3]
            s1 = '.'*start + s1[start:mid1] + '.'*(mid2-mid1) + s1[mid2:end] + '.'*(len(s1)-end)
            s2 = '.'*start + s2[start:mid1] + '.'*(mid2-mid1) + s2[mid2:end] + '.'*(len(s1)-end)
            # print (s1)
            # print (s2)



    start = time.time()
    # Set model details (by Stefan)
    md = RNA.md()
    md.temperature = 37.0
    # md.temperature = 22.0
    # How to treat "dangling end" energies for bases adjacent to helices in free ends and multi-loops.
    # md.dangles = 0
    md.dangles = 2
    md.logML = 0
    md.special_hp = True
    md.noGU = False
    md.noGUclosure = False
    fc = RNA.fold_compound(sequence, md)
    # fc = RNA.fold_compound(sequence)


    bpd = RNA.bp_distance(s1, s2)

    if verbose:
        print(f"Base-pair distance: {bpd}")


    max_energy = None  # Specify upper bound for barrier energy (kcal/mol).
    
    # search_width = None  # Adjust upper bound for findpath search.

    if search_width is None:
        search_width = 2 * bpd

    paths = fc.path_findpath(s1, s2, width=search_width)

    sE = float('-inf')     # saddle energy
    s_pos = 0
    en_min = float('inf')  # minimum, assuming we get below E1


    current_structure = ""

    # compare changes to previous move

    def str_compare(s1, s2):
        for i, (c1, c2) in enumerate(zip(s1, s2)):
            if c1 != c2:
                if c1 == ".":
                    yield i
                else:
                    yield -i

    moves_str = []
    moves_en = []
    moves_pos = []

    for i, step in enumerate(paths):
        moves = [i for i in str_compare('.'+current_structure, '.'+step.s)]

        current_structure = step.s
        moves_str.append(step.s)
        moves_en.append(step.en)

        if moves==[]:
            moves_pos.append((0,0,step.en))
        else:
            moves_pos.append((moves[0],moves[1],step.en))

        if step.en > sE:
            sE = step.en
            s_pos = i


    e1 = round(fc.eval_structure(s1), 2)
    e2 = round(fc.eval_structure(s2), 2)
    barrier = sE - e1

    # minimum between e1 and saddle energy (max energy)
    if len(moves_en[0:s_pos])==0:
        en_min_a = e1
    else:
        en_min_a = min(moves_en[0:s_pos])

    if len(moves_en[s_pos+1:])<=1:
        en_min_b = e2
    else:
        en_min_b = min(moves_en[s_pos+1:])

        # en_min_a = en_min-e1


    end = time.time()
    runtime = end - start


    if verbose:
        print_moves(sequence, s1, s2, moves_pos, move_color='\033[92m')


    return_path = path_class()
    return_path.max_en=sE
    return_path.max_en_pos=s_pos
    return_path.min_en_a=en_min_a
    return_path.min_en_b=en_min_b
    return_path.e1_en=e1
    return_path.e2_en=e2
    return_path.paths = [moves_pos]
    return_path.runtime = runtime

    return return_path
        # return sE, s_pos, b1, b2, e1, e2, moves_str, moves_en


def plot_paths(path, size=500):
    title = f"S: {path.sE:6.2f}  | B: {path.barrier:6.2f} | E[0]:{path.e1:7.2f} E[{len(path.moves_pos)}]:{path.e2:7.2f} | kcal/mol | {path.runtime:6.4f} s"
    # print(path.moves_en)
    plt.plot(path.moves_en)
    plt.ylabel('kcal/mol')
    plt.xlabel('moves')
    plt.title(title)

    return plt

def plot_mpath(*paths, labels=[], outfile='plot.svg', size=500):

    title = ""

    plt.figure()

    for i, path in enumerate(paths):
        data = []
        for p in path:
            data.append(p[2])

        plt.plot(data, label=labels[i])


    plt.ylabel('kcal/mol')
    plt.xlabel('moves')
    plt.title(title)
    plt.legend()

    # plt.show()
    plt.savefig(outfile)

    return plt



if __name__ == '__main__':

    
    sequence = 'UAAAAUGAUCACGGUUUCAGCUUUGGACAGGGCGUUCCACUAAACUCCUGGUGACAUAGAUAUAUUGGAUUGCAACUACUCGUCGGUCCGGUUGGCGUUC'
    s1 = "....(((.((((...............(((((.(((......)))))))))))))))........((((((((........).))))))).........."
    s2 = "..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).))))))))))))).."

    s1 = "....(((.((((...............(((((.(((......)))))))))))))))......(((((((((...........)))))))))........"
    s2 = "..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).))))))))))))).."

    s1 = "....(((.((((...............(((((.(((......))))))))))))))).........(((((((........).))))))..........."
    s2 = "..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).))))))))))))).."

    s1 = ".....((.((((.....(......)..(((((.(((......)))))))))))))).......((((((((((........).)))))))))........"
    s2 = "..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).))))))))))))).."

    s1 = "........((((....((......)).(((((.(((......)))))))))))).........((((((((((........).)))))))))........"
    s2 = "..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).))))))))))))).."

    # s1 = "........((((...............(((((.(((......)))))))))))).........((((((((((........).)))))))))........"
    # s2 = "..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).))))))))))))).."

    # fp call 53 -5.4
    s1 = "..........(((...((((((..((((..((((...(((((......))))).........(((((.....))).))..)))).))))))))))))).."
    s2 = "....(((.((((...............(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........"

    # fp call 49 -3.7
    # s1 = "...........((...((((((..((((..((((...(((((......)))))............((.....))......)))).))))))))))))..."
    # s2 = "....(((.((((...............(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........"

    # fp call 49 -4.9
    # s1 = "..........(((...((((((..((((..((((...((((........))))............((.....))......)))).))))))))))))).."
    # s2 = "....(((.((((...............(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........"

    # fp call 46 -2.29
    s1 = "................(((((..(((((..((((...(((((......)))))..............(....).......)))).))))))))))....."
    s2 = "....(((.((((...............(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........"
    # fp call 46 -2.9
    s1 = "...........((....(((((..((((..((((...(((((......)))))...........................)))).))))))))).))..."
    s2 = "....(((.((((...............(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........"
    # fp call 48 -2.1
    s1 = "..........(((....((((...((((...(((...(((((......)))))...........(((.....))).....)))..)))).)))).))).."
    s2 = "....(((.((((...............(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........"

    # fp call 22 -9.4
    s1 = "..........(((...((((((.....(((((.(((......))))))))................(((((((........).))))))))))))))).."
    s2 = "....(((.((((...............(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........"

    # fp call 6 -15.1
    # s1 = "....(((.((((...((((....))))(((((.(((......)))))))))))))))......(((((((((...........)))))))))........"
    # s2 = "....(((.((((...............(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........"



    sequence = 'UAAAAUGAUCACGGUUUCAGCUUUGGACAGGGCGUUCCACUAAACUCCUGGUGACAUAGAUAUAUUGGAUUGCAACUACUCGUCGGUCCGGUUGGCGUUC'
    s1       = '....(((.((((...............(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........'
    s2       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'
    
    

    # sections
    s1       = '....(((.((((...............(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........'
    s2       = '..........(((...((((((..((((((((.(((......))))))))..............(((.....)))...........))))))))))))..'

    s1       = '..........(((...((((((..((((((((.(((......))))))))..............(((.....)))...........))))))))))))..'
    s2       = '................((((((..((((.(((.(((......))).)))...(((...((((...((.....))..)).))))).)))))))))).....'

    s1       = '................((((((..((((.(((.(((......))).)))...(((...((((...((.....))..)).))))).)))))))))).....'
    s2       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'
    #

    # corrected 7.6 target
    s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'
    s2       = '....(((.(((((((....))).....(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........'
    
    s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'
    s2       = '....(((.((((...((((....))))(((((.(((......)))))))))))))))......((((((((((........).)))))))))........' # mfe? -18.5


    # barriers path
    # s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'
    # s2       = "..........(((...((((((..((((..(((((.(((.........))).)...........................)))).))))))))))))).."

    # s1       = "..........(((...((((((..((((..(((((.(((.........))).)...........................)))).))))))))))))).."
    # s2       = '..........(((...((((((..((((..((((..............((((..((..............))..))))..)))).)))))))))))))..'

    # s1       = '..........(((...((((((..((((..((((..............((((..((..............))..))))..)))).)))))))))))))..'
    # s2       = '..........(((...((((((..((((.(((.(((......))).)))................((.....))...........)))))))))))))..'

    # s1       = '..........(((...((((((..((((.(((.(((......))).)))................((.....))...........)))))))))))))..'
    # s2       = '..........(((...((((((..((((((((.(((......)))))))..((((.(((.....(((.....))))))...)))))))))))))))))..'

    # s1       = '..........(((...((((((..((((((((.(((......)))))))..((((.(((.....(((.....))))))...)))))))))))))))))..'
    # s2       = '..........(((...((((((.....(((((.(((......))))))))................(((((((........).)))))))))))))))..'
    # sequence = "GGAAGCCGGCGAGGCAGUACCAUUAUAUAGUUUGUCUUCCAAGAAUGGGUACGACCGCGGGACCGUUCGGUUAUCGUCUG"
    # s1 = '((.((((((((.((..((((((((....((.....))......))))).)))..)).)).......)))))).)).....'
    # s2 = '((.((((((((.....(((((..........................))))).....)).......)))))).)).....'

    # s1 = '((.((((((((.((..((((((((....((.....))......))))).)))..)).)).......)))))).)).....'
    # s2 = '.(((((((.((.((..(((((..........................)))))..)).)).((...))))))).)).....'

    path = pathfinder(sequence, s1, s2, verbose=True)
    print (path)
    # plt = plot_paths(path)
    # plt.show()

    # print(path)

    # sE, s_pos, b1, b2, e1, e2, moves_str, moves_en = pathfinder(sequence, s1, s2)
    # print(moves_str)
