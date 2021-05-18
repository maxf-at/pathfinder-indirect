#!/usr/bin/env python3
# coding: utf-8

# todo

# implement nb
# gradual nb find at lower search width
# additive approach
# search db when possible %%%
# check barriers

# cyclic dependency direct path moves <-> stabilizing moves

# todo - fix fwd / bwd
# todo - maxE constraint fixed at indirect maximum

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
import pathfinder
import pathfinder_cpp


def print_d(*args):
    """
    debug print for coloredlogs library
    """
    msg = ""
    for m in args:
        msg += str(m)
        msg += " "
    logging.debug(msg)


# https://www.reddit.com/r/Python/comments/27crqg/making_defaultdict_create_defaults_that_are_a/
class key_dependent_dict(defaultdict):
    def __init__(self, f_of_x):
        super().__init__(None)  # base class doesn't get a factory
        self.f_of_x = f_of_x  # save f(x)

    def __missing__(self, key):  # called when a default needed
        ret = self.f_of_x(key)  # calculate default value
        self[key] = ret  # and install it in the dict
        return ret


@dataclass
class Intermediate:
    p_table:      list
    saddle_e:     int
    current_e:    int
    mode:         bool  # 1 = fwd, 0 = bwd
    moves:        list
    energies:     list
    opt:          float
    add_moves:    set

    moves_available: list


def generate_end_s(pt1, moves):

    c_ptable = pt1.copy()

    # generate a random path...

    m = moves.copy()
    # print(m)

    while (len(m)>0):

        
        ls = RNA.loopidx_from_ptable(c_ptable)

        for a, (i,j) in enumerate(m):
            if i > 0 and c_ptable[i] == 0 and c_ptable[j] == 0 and ls[i] == ls[j]:
                c_ptable[i] = j
                c_ptable[j] = i
                break

            if i < 0 and -i == c_ptable[-j]:
                c_ptable[-i] = 0
                c_ptable[-j] = 0
                break
        
        m.remove((i,j))

    return c_ptable



def find_moves(c_ptable, t_ptable, current_add_moves, current_common_moves):
    """
    generator function, yields possible structures 1 move away
    from the original structure by finding fitting i and j with
    RNA pair and loop tables
    c_ptable: current ptable
    end_ptable: s2 end ptable
    """
    # loop table
    ls = RNA.loopidx_from_ptable(c_ptable)

    # extra indirect moves: i and j have to be unpaired, and on the same loop

    # for i, j in self.moves_add:
    for i, j in current_add_moves:
        if i > 0 and c_ptable[i] == 0 and c_ptable[j] == 0 and ls[i] == ls[j]:
            yield i, j, False

        # if i < 0:
        if i < 0 and -i == c_ptable[-j]:
            yield i, j, False


    # print (current_common_moves)

    if len(current_common_moves)>0:

        i = current_common_moves[0][0]
        j = current_common_moves[0][1]
        if i > 0 and c_ptable[i] == 0 and c_ptable[j] == 0 and ls[i] == ls[j]:
            yield i, j, True

        # if i < 0:
        if i < 0 and -i == c_ptable[-j]:
            yield i, j, True

    # for i in range(len(c_ptable)):
    #     if i == 0:
    #         continue

    #     if c_ptable[i] == 0 and t_ptable[i] > i:
    #         j = t_ptable[i]
    #         # found j has to be empty and currently on the same loop as i
    #         if c_ptable[j] == 0 and ls[i] == ls[j]:

    #             # print ("yield", i, j, (i,j) in current_add_moves)
    #             yield i, j
    #     # test for bp removal: i has to be paired with a different j in s2
    #     j = c_ptable[i]
    #     # dont remove things which are present in s2
    #     if c_ptable[i] > i and c_ptable[i] != c_ptable[j] and\
    #             c_ptable[i] != t_ptable[i] and c_ptable[j] != t_ptable[j]:
            
    #         # print ("yield-", -i, -j, (-i,-j) in current_add_moves)
    #         yield -i, -j


def fp_once(sequence, s1, s2, max_energy, width, common_move_list, add_moves=[], mode=True):
    """
    main findpath algorithm (bounded BFS)
    """

    print("322")
    if Debug:
        coloredlogs.DEFAULT_LOG_FORMAT = '%(levelname)s %(message)s'
        coloredlogs.install(level='DEBUG')

    # common_move_list = set(common_move_list)

    # add_moves |= common_move_list


    

    # yield 16 69 False
    # yield- -31 -56 False
    # yield- -32 -55 False

    fc = RNA.fold_compound(sequence)

    pt1 = list(RNA.ptable_from_string(s1))
    pt2 = list(RNA.ptable_from_string(s2))

    e1 = fc.eval_structure_pt(pt1)
    e2 = fc.eval_structure_pt(pt2)

    current_bp = 0
    current_bp_end = RNA.bp_distance(s1, s2)      # current_bp_end = 4
    # paths = [(list(p_tables[s1]), float("-inf"), evals[s1], [])]

    print("bp dist:", current_bp_end, len(add_moves))

    end_p_table = pt2
    init_intermediate = Intermediate(p_table=pt1.copy(), mode=mode, saddle_e=float(
        "-inf"), current_e=e1, moves=[(0, 0, e1)], energies=[0], opt=0, add_moves=add_moves.copy(), moves_available=common_move_list.copy())
    # initial path start with 1 intermediate
    init_path = [init_intermediate]
    # paths start with 1 initial path
    paths = [init_path]

    print_d("start fp once")
    # print_d(add_moves)
    # print_d(common_move_list)

    # dont stop at current_bp_end, consider potential indirect moves
    while (current_bp != current_bp_end):

        # collect all new paths here (next iteration)
        collect_paths = []

        for current_path in paths:
            current_p_table = current_path[-1].p_table
            current_e = current_path[-1].current_e
            current_s = current_path[-1].saddle_e
            # current_string = p_to_s(current_p_table)
            current_moves = current_path[-1].moves
            current_energies = current_path[-1].energies

            current_add_moves = current_path[-1].add_moves.copy()
            current_common_moves = current_path[-1].moves_available.copy()

            # print_d (len(current_add_moves), len(current_common_moves), current_common_moves[0])



            # "try_moves"
            for i, j, flag in find_moves(current_p_table, end_p_table, current_add_moves, current_common_moves):

                # if (i, j) in current_add_moves:
                # continue  # this optional move is already in the path

                # if (i, j) in self.moves_add:
                # current_add_moves.add((i, j))

                # next energy calculations
                next_e = fc.eval_move_pt(
                    current_p_table, i, j) + current_e

                next_e = round(next_e, 2)

                next_p_table = current_p_table.copy()
                if i < 0:
                    next_p_table[-i] = 0
                    next_p_table[-j] = 0
                else:
                    next_p_table[i] = j
                    next_p_table[j] = i

                next_common_moves = current_common_moves.copy()
                if flag:
                    next_common_moves.pop(0)


                # next saddle energy
                next_s = round(max(current_s, next_e), 2)

                # if this move is valid... append to list
                if next_s <= max_energy:

                    next_moves = current_moves.copy()
                    next_moves.append((i, j, next_e))

                    next_energies = current_energies.copy()
                    if next_e < e1:
                        next_energies.append(round(next_e-e1, 2))
                    else:
                        next_energies.append(0)

                    new_intermediate = Intermediate(p_table=next_p_table, mode=mode, saddle_e=next_s, current_e=next_e,
                                                    moves=next_moves, energies=next_energies, opt=[], add_moves=current_add_moves, moves_available=next_common_moves)

                    new_path = current_path.copy() + [new_intermediate]
                    collect_paths.append(new_path)

        # first sorting step
        collect_paths.sort(key=lambda x: (x[-1].p_table, x[-1].saddle_e))

        last_ptable = []
        # last_ptable = collect_paths[-1][0].p_table
        # last_ptable = collect_paths[0][-1].p_table

        print_d(current_bp, "sort done")

        # remove duplicates ptables
        if current_bp+1 != current_bp_end:
            for i in range(len(collect_paths)):
                c_ptable = collect_paths[i][-1].p_table
                if c_ptable == last_ptable:
                    # set saddle energy high
                    collect_paths[i][-1].saddle_e = 999
                else:
                    last_ptable = c_ptable
            i = 0
            while True:
                if i > len(collect_paths)-1:
                    break
                if collect_paths[i][-1].saddle_e == 999:
                    collect_paths.pop(i)
                    continue
                i += 1

        # second sorting step
        collect_paths.sort(key=lambda x: (x[-1].saddle_e, x[-1].current_e))

        # print("iteration", current_bp, len(collect_paths))

        # discard paths
        collect_paths = collect_paths[:width]

        for p in collect_paths[:5]:
            a = p[-1].saddle_e
            cs = RNA.db_from_ptable(p[-1].p_table)

            print(cs, a)

        # return valid paths if we're past the minimum bp_dist
        if current_bp >= current_bp_end-1:  # or current_bp==30:
            for i in range(len(collect_paths)):
                if collect_paths[i][-1].p_table == list(end_p_table):
                    # print("y")
                    yield collect_paths[i][-1]
                # yield collect_paths[i][-1]

        # next iteration
        paths = collect_paths
        current_bp += 1

    # return remaining paths
    if paths:
        for path in paths:
            if path:
                yield path[-1]


def find_path(sequence, s_ref, destination_list, max_energy=float("inf"), results=1, search_width=1000, Debug=False, Verbose=False):
    """
    indirect findpath, main function

    settings:
    indirect_iterations, default value 2 (means 1 direct pass, 1 indirect pass)

    """

    if Debug:
        coloredlogs.DEFAULT_LOG_FORMAT = '%(levelname)s %(message)s'
        coloredlogs.install(level='DEBUG')

    s1 = s_ref
    pt1 = list(RNA.ptable_from_string(s1))
    # generate move list

    destination_list.sort()

    moves_dest = [set() for i in destination_list]
    moves_dest_add = [set() for i in destination_list]

    bp_dist_dest = [None for i in destination_list]
    pt_dest = [None for i in destination_list]

    min_bp_dist = float("inf")
    max_bp_dist = float("-inf")

    print_d(s_ref)
    print_d("destinations:")

    for a, s in enumerate(destination_list):
        bp_dist = RNA.bp_distance(s, s_ref)

        current_move_list = set()

        c_pt = list(RNA.ptable_from_string(s))
        pt_dest[a] = c_pt

        for i in range(1, c_pt[0]):
            if pt1[i] != c_pt[i]:
                if i < pt1[i]:
                    # bp deletion
                    current_move_list.add((-i, -pt1[i]))

                if i < c_pt[i]:
                    # add bp
                    current_move_list.add((i, c_pt[i]))

        moves_dest[a] = current_move_list.copy()

        bp_dist = len(current_move_list)
        bp_dist_dest[a] = bp_dist

        if bp_dist < min_bp_dist:
            min_bp_dist = bp_dist

        if bp_dist > max_bp_dist:
            max_bp_dist = bp_dist
        print_d("dest:", a, s, bp_dist)
        print_d(len(current_move_list), current_move_list)
        # break

    # One-Liner to intersect a list of sets

    all_moves = moves_dest[0].union(*moves_dest)
    common_moves = moves_dest[0].intersection(*moves_dest)

    print_d("common moves:", common_moves, len(common_moves))

    additional_moves = set()

    for a in range(len(moves_dest)):
        for current_set in moves_dest[a]:

            if current_set not in common_moves:

                # add inverse move
                i = -current_set[0]
                j = -current_set[1]
                moves_dest_add[a].add((-i, -j))
                # additional_moves.add((-a,-b))

        print_d("non-common moves", a, moves_dest_add[a], len(moves_dest_add[a]))

    # generate s_common from intersection moves

    pt_common = generate_end_s(pt1, common_moves)


    s_common = RNA.db_from_ptable(pt_common)
    print ("bp dist common", RNA.bp_distance(s1, s_common), len(common_moves))

    print_d(s1)
    print_d(s_common)

    search_width = 20
    # Verbose = False
    Verbose = True

    r = pathfinder_cpp.find_path(sequence, s1, s_common, search_width=search_width, Verbose=Verbose)

    common_move_list = list(reversed(r.paths))[0:-1]  # back to front, skip the 0,0 move
    common_move_list = r.paths[1:]


    print("common moves:", len(common_move_list), common_move_list)

    add_moves = moves_dest_add[0]

    print("add moves", len(add_moves), add_moves)

    print (moves_dest[a]-common_moves)

    pt_d0 = generate_end_s(pt_dest[0], common_moves)
    s_d0 = RNA.db_from_ptable(pt_d0)

    s0 = destination_list[0]

    r0 = pathfinder_cpp.find_path(sequence, s_common, s0, search_width=search_width, Verbose=Verbose)

    # r = pathfinder_cpp.find_path(sequence, s1, s0, search_width=search_width, Verbose=Verbose)
    # return
    # import merge_cpp
    # merge_search_width = 10
    # merge_cpp.merge_findpath(sequence, s1, s0, r, r0, [],[], Debug=Debug, Verbose=Verbose, merge_search_width=merge_search_width)

    print("fp once")

    all_paths = []
    for path in fp_once(sequence, s1, s0, max_energy, search_width, common_move_list, add_moves, True):
        all_paths.append(path)

    all_paths.sort(key=lambda x: x.saddle_e)

    # print (all_paths[0])

    print_moves(sequence, s1, s2, all_paths[0].moves)


    return

    intersection = intersection.union(additional_moves)
    # print(intersection)
    # print ("add moves:", additional_moves)

    fp_class.moves_add = intersection.copy()
    # add inverse moves to undo bs...

    for k, (s, bp_dist) in enumerate(fp_class.destinations):
        for a, b in additional_moves:
            if (-a, -b) in fp_class.available[k]:
                fp_class.available[k].add((a, b))
                # print (-a,-b)

    # destinations: start with the largest bp dist
    # print (fp_class.destinations)

    # all_moves[0] = move_list # combined

    # fp_class.available = np.zeros([len(destination_list), len(move_list)]).astype(bool)
    # # which moves are available where
    # for b, t in enumerate(move_list):
    #     for a in range(len(destination_list)):
    #         fp_class.available[a][b] = t in all_moves[a]

    # fp_class.available = all_moves
    # print (fp_class.available)

    # max_en = 9999999999

    print_d("~FWD~~~~~~")
    result_dict = defaultdict(lambda: 9999)

    # fwd
    path_generator = fp_class.find_path_once(
        start_list, destination_list, max_energy, current_width, mode=True, sort_min=False, Debug=Debug, Verbose=Verbose)

    for path in path_generator:
        max_en = path.saddle_e

        end_s = RNA.db_from_ptable(path.p_table)  # ?

        # if end_s not in destination_list:
        # continue

        start_index = path.start_pt
        start_s = fp_class.destinations[start_index][0]

        # print_moves(sequence, s1, start_s, path.moves)

        if max_en < result_dict[end_s]:
            result_dict[end_s] = max_en

    for key in result_dict:
        print("fwd:", key, result_dict[key])
# #

    if Debug:
        coloredlogs.DEFAULT_LOG_FORMAT = '%(levelname)s %(message)s'
        coloredlogs.install(level='DEBUG')

    print_d("~BWD~~~~~~")
    max_energy = 999
    # result_dict = defaultdict(lambda: 9999)

    # path_generator = fp_class.find_path_once(
    #     start_list, destination_list, max_energy, current_width, mode=False, sort_min=False, Debug=Debug, Verbose=Verbose)

    # for path in path_generator:
    #     max_en = path.saddle_e

    #     end_s = RNA.db_from_ptable(path.p_table)

    #     start_index = path.start_pt
    #     start_s = fp_class.destinations[start_index][0]

    #     # start_s = 3
    #     # print (start_s)

    #     # print ("f", start_s, max_en, path.moves)
    #     # print_moves(sequence, start_s, end_s, path.moves)

    #     for e in path.embedded:
    #         embedded_en = path.embedded[e]
    #         if embedded_en <= max_en:
    #             em_s = fp_class.destinations[e][0]
    #             if max_en < result_dict[em_s]:
    #                 result_dict[em_s] = max_en
    #             # print ("e", e, path.embedded[e])

    #     if max_en < result_dict[start_s]:
    # result_dict[start_s] = max_en

    # BWD END
    return_energy = float("inf")

    for key in result_dict:
        print("final:", key, result_dict[key])
        if result_dict[key] < return_energy:
            return_energy = result_dict[key]

    # code here to not just return nothing
    return return_energy


if __name__ == '__main__':

    sequence = 'CCGCGCCAGCGACGGAUUGCUCACUAAAUCUACAGGAGUUAUAGAUCAGCAUCCGGACUAAUGAGACGUGGCAUUUCCACUGGUUAUAUCACUGGCGAGUGCCAAUAACUACCUUCGUAGUACGGUAGAGGCCCUUCUACAUCCACGCCU'
    s1 = '.(((....))).(((((.(((..(((..(((...)))....)))...)))))))).............(((((((((((.((((...)))).))).))))))))...(((((....)))))...((((((....))))))..........'
    s2 = '.(((....))).((((.((((..(((((.((.....)))).)))...))))))))........((...(((((((((((.((((...)))).))).))))))))...(((((....)))))...((((((....))))))........))'
    intermediates = ['.(((....))).((((.((((..(((..(((...)))....)))...))))))))...........(((((((((((((.((((...)))).))).)))))))....(((((....))))))))((((((....))))))..........', '.(((....))).((((.((((..(((...((.....))...)))...))))))))...........(((((((((((((.((((...)))).))).)))))))....(((((....))))))))((((((....))))))..........',
                     '.(((....))).((((.((((..(((.....((....))..)))...)))))))).............(((((((((((.((((...)))).))).))))))))...(((((....)))))...((((((....))))))..........']

    # 15
    sequence = 'AAAACGGCUAGACGUUUAUGACUGGCGAUAAUUAUUUGCAUGAUGGAUCCAAUCCAAUGCAGUCGAGUGAAGUACCGUCGUAGUAGCUAUAACCACCUUAUUCUUGCCCCCUGCAAGCGCUAGAUACUGAUCCAUGCUAAUAAAGUAGUC'
    s1 = '.....(((((.....(((((((.((.....(((..((((((..(((((...)))))))))))...)))......))))))))))))))...((.((.((((((((((.....))))).((..(((....)))...)).))))).)).)).'
    s2 = '.(((((......)))))....((((((........((((((..(((((...)))))))))))....(((.(((..(......)..))).....)))......(((((.....))))))))))).....(((...((((.....)))))))'
    intermediates = ['.(((((......)))))((..((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).)))))))).........((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).))))))..(((.......(((.....)))))).', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).))))))....((...))((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).)))))).....(((...((((.....)))))))',
                     '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).))))))...........((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..(((((.((....((.((((...)))).))...)))))))((((.....)))).))))))...........((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((.(((....((.((((...)))).))...)))))))((((.....)))).))))))...........((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))...(((((.((...((.((((...)))).)).)).)))))(((((.....)))))))))))...........((((.....))))...', '.(((((......))))).....(((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).)))))(((....)))..((((.....))))...']
    max_energy = -12.9
    max_energy = -14.0

    # 150 nt, -20.10 kcal/mol direct, indirect -24.50
    # sequence = 'AAAAUAAUGUACCGGACAUUCGCGCACGACCACCAUAUGGCAGAGCAUGUGUCUGUGGACCCACUAUAGCUGGGGCGCUUAACCCCAGAAAAGUAUCUUCGGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAA'
    # s1       = '.............((((((..((...(..(((.....)))..).))..))))))(((((((....(((.((((((.......)))))).....)))....)))))))(((......((.((.........)).))......)))......'
    # s2       = '......(((.((((((...(((....))).(((((((((......))))))...)))((...(((....((((((.......))))))...))).)))))))).)))(((......((.((.........)).))......)))......'
    # intermediates = ['......(((.((((((.(((((....)))((((((((((......))))))...))))....(((....((((((.......))))))...))))).)))))).)))(((......((.((.........)).))......)))......', '......(((.((((((...(((....)))((((((((((......))))))...))))....(((....((((((.......))))))...)))...)))))).)))(((......((.((.........)).))......)))......', '......(((.((((((.............((((((((((......))))))...))))....(((....((((((.......))))))...)))...)))))).)))(((......((.((.........)).))......)))......', '..........((((((...(((....)))((((((((((......))))))...))))....(((....((((((.......))))))...)))...))))))....(((......((.((.........)).))......)))......']

    # sequence = 'AUCAUUGCCUAAGGGACAUUUCCACGUGUAAGAAACCCGGAAACAUCGUUUGGCUUAGAAUGAACAGAUAGAUUGCCUACAACAAAAUACUGAUGACGGGGGAGCCGUGACAUGGCGGCCACCGCUAGUCGGCCUACGCUUUAAAACAUA'
    # s1       = '.((((.(((.(((.((..(((((..((.......))..)))))..)).))))))..........(((.((..((..........)).)))))))))((..((.((((.....((((((...))))))..)))))).))............'
    # s2       = '.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).(((.....(((.......)))....)))..((((.((..(((((......)))))..)).)..)))(((....)))..........'

    # intermediates = [
    # '.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    # '.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).........(((....))).......(((((..(((((..(((((...)))))..)).)))...)))))..................',
    # '.((((((((.(((.((..(((((..(.((.....))).)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    # '.((((((((.(((.((..(((((...............)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    # '.((((((((.(((.((..(((((...............)))))..)).))))))....)))))..........................(((((..(((((..(((((...)))))..)).)))...)))))..................',
    # '......(((.(((.((..(((((..((.......))..)))))..)).))))))..............(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................'
    # ]

    # simple 60 nt ex
    sequence = "GAAAGUUUUCAAGAAAGUAAUAUAGUGUCCCAGUCAGGGGGUAGUGGCGCGUUCUGGACGCAUAGUAUUC"
    s1 = "(((....))).........(((((((((.((((...(.(.((....)).).).)))))))).)).))).."
    s2 = ".................((.(((.((((.((((..(.(..((....)).).).))))))))))).))..."
    intermediates = ['(((....))).....((((.(((.(((((.....(((((.((......)).))))))))))))).)))).', '...............((((.(((.((((((.....((((.((......)).))))))))))))).)))).', '...............((((.(((.(((((.((....(((.((......)).))))))))))))).)))).', '....................(((.((((((.....((((.((......)).)))))))))))))......']

    # todo: bp_min from min_bp_dist onwards - best candidate
    # bwds + indirect

    # search_width = 4
    # search_width = 20
    search_width = 150
    # search_width = 500

    Verbose = True
    Debug = True
    # Debug = False

    # max_energy = -13.5
    # max_energy = 999

    intermediates = list(set(intermediates))

    en = find_path(sequence, s1, intermediates, max_energy=max_energy,
                   search_width=search_width, Debug=Debug, Verbose=Verbose)
    print("to s1", en)

    # en = find_path(sequence, s2, intermediates, max_energy=max_energy,
    #                   search_width=search_width, Debug=Debug, Verbose=Verbose)
    # print("to s2", en)

    search_width = 500
    search_width = 20

    print('orig findpath:')

    Verbose = False
    # Verbose = True

    destinations1 = [[i, 0] for i in intermediates]
    for a, s in enumerate(intermediates):
        bp_dist = RNA.bp_distance(s, s1)
        destinations1[a] = (s, bp_dist)
    destinations1.sort(key=lambda x: -x[1])

    destinations2 = [[i, 0] for i in intermediates]
    for a, s in enumerate(intermediates):
        bp_dist = RNA.bp_distance(s, s2)
        destinations2[a] = (s, bp_dist)
    destinations2.sort(key=lambda x: -x[1])

    for number, (i, _) in enumerate(destinations1.copy()):
        print("s1: ~~~~~~~~~~")
        # print (number, i, _)
        r1 = pathfinder.pathfinder(
            sequence, s1, i, search_width=search_width, verbose=Verbose)
        print(number, i, r1.max_en, RNA.bp_distance(s1, i))

    # for number, (i,_) in enumerate(destinations2.copy()):
    #     print ("S2: ~~~~~~~~~~")
    #     # print (number, i, _)
    #     r2 = pathfinder.pathfinder(
    #         sequence, s2, i, search_width=search_width, verbose=Verbose)
    #     print(number, i, r2.max_en, RNA.bp_distance(s2,i))

    # pathfinder.pathfinder(sequence, s1, s2, search_width=search_width, verbose=True)


def detect_local_minimum(fc, structure):
    # perform gradient walk from sample to determine direct local minimum
    pt = RNA.IntVector(RNA.ptable(structure))
    fc.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
    return RNA.db_from_ptable(list(pt))
