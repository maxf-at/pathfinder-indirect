#!/usr/bin/env python3
# coding: utf-8

import sys
import random
import string
import RNA
import numpy as np
import os
import subprocess
from dataclasses import dataclass, field

import matplotlib.pyplot as plt
from PIL import Image
from IPython.display import SVG, display


@dataclass
class path_class:
    min_en_a: float = float("inf")
    min_en_b: float = float("inf")
    max_en: float = float("-inf")
    max_en_pos: int = 0
    

    e1_en: float = 0
    e2_en: float = 0

    dist: int = 0

    s1: str = ""
    s2: str = ""

    paths: list = field(default_factory=list)
    section: list = field(default_factory=list)
    G: list = field(default_factory=list)

    runtime: float = 0
    runtime_ext: float = 0
    runtime_merge: float = 0
    runtime_merge_ext: float = 0


def p_to_s(tb):
    """ptable --> string"""
    s = ""
    for i, p in enumerate(tb):
        if i == 0:
            continue
        if p == 0:
            s += "."
        elif p > i:
            s += "("
        else:
            s += ")"
    return s


"""
read and save functions to save a sequence name,
structure, s1, s2 into a .seq file
"""


def save_seq(sequence, s1, s2, filename):
    if len(sequence) != len(s1) and len(s1) != (s2):
        raise
    with open(filename, 'w') as f:
        f.write(f'>{filename}\n')
        f.write(f'{sequence}\n')
        f.write(f'{s1}\n')
        f.write(f'{s2}\n')


def read_seq(filename):
    # check if exists etc...
    with open(filename, 'r') as f:
        return_list = []
        for line in f:
            line = line.strip()
            if line[0] == ">":
                continue
            return_list.append(line)
        return return_list[0], return_list[1], return_list[2]


def generate_structures(length=150):

    subopt_list = []
    while len(subopt_list) < 3:
        sequence = RNA.random_string(length, 'AUGC')
        cmd = f'printf "{sequence}" | RNAsubopt -p 2'
        result = subprocess.check_output(cmd, shell=True, encoding="utf8")
        subopt_list = result.split()
        break

    sequence, s1, s2 = subopt_list

    fc = RNA.fold_compound(sequence)

    # pt1 = RNA.ptable_from_string(s1)
    # pt2 = RNA.ptable_from_string(s2)
    # print (pt1)
    # test = fc.path_gradient(fc, (1,pt1))

    return sequence, s1, s2


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

    fc = RNA.fold_compound(sequence)
    e1 = en = round(fc.eval_structure(s), 2)
    max_en = float("-inf")



    output_rows = []
    moves_i_j = [(x[0],x[1]) for x in moves]

    # preprocessing - generate strings & energies if required
    for a, move in enumerate(moves):
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
            colored_s = s[0:pos_i] + c.RED + c.BOLD + s[pos_i] + c.ENDC +\
                s[pos_i+1:pos_j] + c.RED + c.BOLD + \
                s[pos_j] + c.ENDC + s[pos_j+1:]
            info = f'{c.RED}[{i:4}, {j:4} ]{c.ENDC}'
        else:  # direct move
            colored_s = s[0:pos_i] + move_color + c.BOLD + s[pos_i] + c.ENDC +\
                s[pos_i+1:pos_j] + move_color + c.BOLD + \
                s[pos_j] + c.ENDC + s[pos_j+1:]
            info = f'{move_color}[{i:4}, {j:4} ]{c.ENDC}'

        # mark pos x
        # x = 10
        # colored_s = colored_s[0:x] + c.CYAN + colored_s[x:x+1] + c.ENDC + colored_s[x+1:]

        if en == max_en:
            info += f' {c.RED}{c.BOLD}{en:6.2f}{c.ENDC}'
        else:
            info += f' {en:6.2f}'

        if Verbose:
            if include != None:
                if abs(i) in include:
                    print(f"{info}")
            elif exclude != None:
                if abs(i) not in exclude:
                    print(f"{info}")
            else:
                print(f"{colored_s} {info}")

    barrier = max_en - e1
    if Verbose: print(
        f"S: {max_en:6.2f} kcal/mol | B: {barrier:6.2f} kcal/mol | E[start]:{e1:6.2f} E[end]:{e2:6.2f}")
    
    return max_en


def plot_layout(sequence, structure, size=500, layout=0):
    # echo -e
    cmd = f'printf ">test\n{sequence}\n{structure}" | RNAplot -t {layout}'
    # print(cmd)
    subprocess.check_output(cmd, shell=True, encoding="utf8")
    img = Image.open('test_ss.ps').resize((size, size), Image.ANTIALIAS)
    # img = Image.open('test_ss.ps')
    # os.remove("test_ss.ps")
    return img


def combine_plots(img_list, size):
    new_img = Image.new('RGB', (size*len(img_list), size))
    offset = 0
    for img in img_list:
        new_img.paste(img, (offset, 0))
        offset += size
    return new_img


def plot_2_seq(sequence, s1, s2, size=500, layout=0):
    s1_plot = plot_layout(sequence, s1, layout=layout, size=size)
    os.rename(r'test_ss.ps', r'plot_1.ps')
    s2_plot = plot_layout(sequence, s2, layout=layout, size=size)
    os.rename(r'test_ss.ps', r'plot_2.ps')
    return combine_plots([s1_plot, s2_plot], size=size)


def print_spacer(s):

    for i in range(len(s)):
        if i % 10 == 0:
            print("|", end='')
        else:
            print(".", end='')
    print()


def print_tables(s, s2, tb1, tb2):
    """
    function to print pair tables, loop tables, etc.
    """

    # https://stackoverflow.com/questions/13673060/split-string-into-strings-by-length
    def split2len(s, n):
        def _f(s, n):
            while s:
                yield s[:n]
                s = s[n:]
        return list(_f(s, n))

    s = "?" + s
    seqs1 = split2len(s, 40)
    s2 = "?" + s2
    seqs2 = split2len(s2, 40)
    tbs1 = split2len(tb1, 40)
    tbs2 = split2len(tb2, 40)

    # offset = 0

    for i in range(len(seqs1)):

        offset = 40 * i
        s1 = seqs1[i]
        s2 = seqs2[i]
        tb1 = tbs1[i]
        tb2 = tbs2[i]

        # for arg in args:
        # 	if len(s) != len(arg):
        # 		return

        for i in range(len(s1)):
            print(f'{i+offset:3}', end=" ")
        print()
        # print()

        print(" ", end="")
        for char in s1:
            print(f'{char:3}', end=" ")
        print()

        print(" ", end="")
        for char in s2:
            print(f'{char:3}', end=" ")
        print()

        # print()

        for char in tb1:
            # if char != 0: char -= 1
            print(f'{char:3}', end=" ")
        print()
        for char in tb2:
            # if char != 0: char -= 1
            print(f'{char:3}', end=" ")
        print()

        print("-------------------")


# build a mountain plot from 2 structures

def mountain_plot(sequence, s1, s2):

    ptables_s1 = RNA.ptable_from_string(s1)
    ptables_s2 = RNA.ptable_from_string(s2)

    ltables_s1 = RNA.loopidx_from_ptable(ptables_s1)
    ltables_s2 = RNA.loopidx_from_ptable(ptables_s2)

    compatible = [None]*len(ptables_s2)

    for i in range(len(ptables_s1)):
        if i == 0:
            continue
        if ptables_s1[i] == ptables_s2[i]:
            compatible[i] = -3

    fig = plt.figure(figsize=(12, 7), dpi=80)
    ax = plt.axes()

    plt.plot(ltables_s1[1:], color='blue')
    plt.plot(ltables_s2[1:], color='orange')
    plt.plot(compatible, color='red')

    return plt


def remove_constant_sections(sequence, s1, s2):

    # if s1[0] == s2[0]:
    #     compatible = [1] + [0]*(len(sequence)-1)
    # else:
    #     compatible = [0] + [0]*(len(sequence)-1) # 0 = "."

    remove = [0] * len(sequence)  # 0 = "."

    fc = RNA.fold_compound(sequence)
    p1 = RNA.ptable_from_string(s1)
    p2 = RNA.ptable_from_string(s2)

    # for iter, (i,j) in enumerate(zip(p1, p2)):
    #     if i==j and i!=0:
    #         compatible[iter] = 1

    def last_bp(i):
        start_i = i
        i -= 1
        while (p1[i] == 0 and p2[i] == 0 and i >= 0):
            i -= 1
        # print ("start j at", start_i, p1[start_i], "to", i,p1[i]  )
        if p1[i] == p2[i] and i > p1[i]:
            return i
        return False
        # return i

    def next_bp(i):
        start_i = i
        i += 1
        while (p1[i] == 0 and p2[i] == 0 and i+1 < len(sequence)):
            i += 1
        # print ("start i at", start_i, p1[start_i], "to", i,p1[i]  )
        if p1[i] == p2[i] and i < p1[i]:
            return i
        return False

    # if i==( and i+1==( --> ...

    seq_trim = ""
    s1_trim = ""
    s2_trim = ""

    last_i = None
    last_j = None

    for i in range(len(sequence)):

        i += 1
        # if p1[i]==0 and p2[i]==0 and last_bp(i) and next_bp(i):
        #     compatible[i] = 1

        # remove closing brackets if compat.
        if (p1[i] == 0 and p2[i] == 0) or (p1[i] == p2[i] and i > p1[i]):
            continue
        if p1[i] == p2[i] and i < p1[i]:  # opening brackets only
            j = p1[i]
            if not last_i:
                last_i = i
                last_j = j
            # lookahead: see if the next i / next j is incompatible,
            # if the next i or j is incompatible, mark the section last_i:i

        if last_i and last_j and not next_bp(i) and not last_bp(j):
            # print ("here", last_i, next_bp(last_i), "->", i, j, "<-", last_bp(last_j), last_j)

            if next_bp(last_i) == i or j == last_bp(last_j):
                # print ("ignore")
                last_i = None
                last_j = None
                continue

            # remove[next_bp(last_i):i] = [1] * (i-next_bp(last_i))
            # remove[j:last_bp(last_j)] = [1] * (last_bp(last_j)-j)

            remove[last_i-1:i-1] = [1] * (i-last_i)
            remove[j-1:last_j-1] = [1] * (last_j-j)
            # print ("".join([str(elem) for elem in remove]))

            last_i = None
            last_j = None

        if p1[i] != p2[i]:
            last_i = None
            last_j = None

    for i in range(len(sequence)):
        if remove[i] == 0:
            seq_trim += sequence[i]
            s1_trim += s1[i]
            s2_trim += s2[i]

    # print (sequence)
    # print (s1)
    # print (s2)
    # print ("".join([str(elem) for elem in remove]))
    # print (seq_trim)
    # print (s1_trim)
    # print (s2_trim)

    # fc_t = RNA.fold_compound(seq_trim)
    # delta_G = fc.eval_structure(s1) - fc_t.eval_structure(s1_trim)
    # return seq_trim, s1_trim, s2_trim, delta_G

    return seq_trim, s1_trim, s2_trim
    # print (sequence[compatible])


if __name__ == '__main__':
    # test = generate_structures()
    # for s in test:
    #     print (s)

    sequence = "AAAUGAUUCAACUUUUCAAAAAAUUAAUAACAUUUUCUCUUAUCGAGAGUUGGGCGAGGGAUUGGCCUUUUGACCCCAACAGCAACCGACCGUAAUACCAUUGUGAAAUGGGGCGCACUGCUUUUCGCGCCGAGACUGAUGUCUCAUAAGGCACGGUGCUAAUUCCAUCAGAUUGUGUCUGAGAGAUGAGAGAGGCAGUGUUUUACGUAGAAAAGCCUCUUUCUCUCAUGGGAAAGAGGCUUUUUGU"
    s1 = "............................((((((((((((((((....((((((((((((......))))))..)))))).(((.(((.((......(((((....)))))(((((..........)))))(((((....)))))....))..))))))........((((((...))))))..))))))))))..)))))).......((((((((((((((((....)))))))))))))))).."
    s2 = "............................(((((((((((((((((((.((((((.(((((.....)))))....))))))(((.((((.((......(((((....)))))(((((..........)))))(((((....)))))....))..)))))))..)))..((((((...))))))..)))))))))))..))))).......((((((((((((((((....)))))))))))))))).."

    # sequence = "ACACANUCAGAUUUCCUNNGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUAAGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUU"
    # s1 =       "..(((.((((.....))..)))))...(((((((((((((.......))))))).)))))).((((((((((.....))))))))))."
    # s2 =       "..(((.((((.....))..)))))................(...(((((....)))))...)((((((((((.....))))))))))."

    sequence = "UGCUAGUCAGGUUGUCAGUCUGUUUUUCGACCCGGACCCCGUAAAGACCCUGGCAAACAUUUCCCCGAAACACGUCGCGGGACGUAAGUCUGAGGCGCCAUGCGAAAACCAGACCGCUGGAGCGACGACCUCAGGAAACUAUGGAGACUU"
    s1 = "((((((...((((..(.(((((..........)))))...)....))))))))))......(((((((......))).))))((((..(((((((((...(((.....((((....)))).))).)).)))))))....))))......."
    s2 = "....((((.(..((((((...((((((.....(((...))).)))))).))))))..)...(((((((......))).))))((((.((((((((((((..(((..........)))..).)))....))))))...))))))..))))."

    print_tables(s1, s2, RNA.ptable_from_string(
        s1), RNA.ptable_from_string(s2))
    print(remove_constant_sections(sequence, s1, s2))
