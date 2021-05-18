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

    embedded: dict
    start_pt: list



class Fp_class():
    def __init__(self, sequence, s1, s2, add_moves=[]):

        self.fc = RNA.fold_compound(sequence)
        self.sequence = sequence
        self.s1 = s1
        self.s2 = s2
        self.destinations = dict()
        self.min_bp_dist = 0
        self.max_bp_dist = 0

        # energy caching
        self.evals = key_dependent_dict(
            lambda x: round(self.fc.eval_structure(x), 2))
        self.p_tables = key_dependent_dict(
            lambda x: RNA.ptable_from_string(x))
        self.l_tables = lambda x: RNA.loopidx_from_ptable(x)
        self.bp_dist = lambda x, y: RNA.bp_distance(x, y)

        # additional class variables for indirect paths
        self.processed_stacks = set()
        self.moves_add = add_moves
        self.available = dict()

    def find_moves(self, c_ptable, current_sets, current_start_pt, mode):
        """
        generator function, yields possible structures 1 move away
        from the original structure by finding fitting i and j with
        RNA pair and loop tables
        c_ptable: current ptable
        end_ptable: s2 end ptable
        """
        # loop table
        ls = self.l_tables(c_ptable)
        
        if mode:
            # for a in current_sets:        
            #     cs = self.available[a]
            #     # print ("cs", cs)
            #     # for i, j in self.moves_add:
            #     for i, j in cs:
            #         if i > 0 and c_ptable[i] == 0 and c_ptable[j] == 0 and ls[i] == ls[j]:
            #             yield i, j

            #         # if i < 0:
            #         if i < 0 and -i == c_ptable[-j]:
            #             yield i, j
            for i, j in self.moves_add:
                if i > 0 and c_ptable[i] == 0 and c_ptable[j] == 0 and ls[i] == ls[j]:
                    yield i, j

                # if i < 0:
                if i < 0 and -i == c_ptable[-j]:
                    yield i, j
        

        if not mode:
            for i, j in self.available[current_start_pt]:

                i = -i
                j = -j
                if i > 0 and c_ptable[i] == 0 and c_ptable[j] == 0 and ls[i] == ls[j]:
                    yield i, j
                # if i < 0:
                if i < 0 and -i == c_ptable[-j]:
                    yield i, j


    def find_path_once(self, start_s, destination_list, max_energy, width, mode=True, sort_min=False, Debug=False, Verbose=False):
        """
        main findpath algorithm (bounded BFS)
        """



        # caching
        evals = self.evals
        p_tables = self.p_tables

        # add_sequences = [[]] * self.max_bp_dist*2
        add_sequences = defaultdict(list)

        s1 = ""
        s2 = ""

        if (mode):  # 1 start, many dests
            # print ("start fwd")
            s1 = start_s
            s2 = destination_list[0]

            # print (s1)
            # print (s2)

        else:
            # pass
            # print ("start bwd")

            # sort by bp dist
            # a = sorted(self.destinations.items(), key=lambda x: -x[1])

            s1 = self.destinations[0][0] # s1 has the largest bp_dist
            s2 = start_s

            # print ("s1 bwd:", s1, 0)
            number = 1
            for seq, c_bp in self.destinations[1:]:
                bp_offset = self.max_bp_dist-c_bp
                # print ("de bwd:", number, seq, bp_offset)
                add_sequences[bp_offset].append((seq,number))
                number += 1


            # print (add_sequences[0])
            # return

        current_bp = 0
        # current_bp_end = bp_dist(s1, s2)      # current_bp_end = 4

        s1_en = round(evals[s1], 2)

        init_intermediate = Intermediate(p_table=list(p_tables[s1]), mode=mode, saddle_e=float(
            "-inf"), current_e=evals[s1], moves=[(0, 0, s1_en)], energies=[0], opt=0, add_moves=set(), start_pt=0, embedded={})
        # initial path start with 1 intermediate
        init_path = [init_intermediate]
        # paths start with 1 initial path
        paths = [init_path]

        # add_sequences[1] = ["..........(((...((((((..((((..((((...((((........))))...........(((.....))).....)))).))))))))))))).."]

        # dont stop at current_bp_end, consider potential indirect moves
        while (current_bp != self.max_bp_dist+20):

            if not mode and add_sequences[current_bp]:
                for id_s, (add_seq,number) in enumerate(add_sequences[current_bp]):

                    available = False
                    #check if the added sequence is already in paths
                    for path_index, current_path in enumerate(paths):
                        current_p_table = current_path[-1].p_table
                        current_string = p_to_s(current_p_table)
                        if current_string == add_seq:
                            available = path_index
                            break

                    if available:
                        current_e=evals[add_seq]                        
                        prev_s = paths[path_index][-1].saddle_e                                                
                        
                        if prev_s > current_e:
                            # change around---                            
                            paths[path_index][-1].embedded[paths[path_index][-1].start_pt] = prev_s
                            paths[path_index][-1].start_pt = number
                            paths[path_index][-1].saddle_e = current_e
                            paths[path_index][-1].current_en = current_e
                        else:
                            # add embedding                            
                            paths[path_index][-1].embedded[number] = current_e

                        # paths[path_index][-1].embedded[number] = current_e
                        
                        print_d ("add em", number, add_seq, prev_s, current_e)
                        # paths[path_index][-1].embedded[number] = -999
                        continue

                    print_d ("add n ", number, add_seq)
                    init_intermediate = Intermediate(p_table=list(p_tables[add_seq]), mode=mode, saddle_e=evals[add_seq], current_e=evals[add_seq], moves=[(0, 0, s1_en)], energies=[0], opt=0, add_moves=set(), start_pt=number, embedded={})
                    init_path = [init_intermediate]
                    paths.append(init_path)

            # collect all new paths here (next iteration)
            collect_paths = []

            for current_path in paths:
                current_p_table = current_path[-1].p_table
                current_e = current_path[-1].current_e
                current_s = current_path[-1].saddle_e
                current_string = p_to_s(current_p_table)
                current_moves = current_path[-1].moves
                current_energies = current_path[-1].energies

                # debug_s = RNA.db_from_ptable(current_p_table)[0:85]
                # print_d(debug_s, current_bp)

                current_start_pt = current_path[-1].start_pt
                current_add_moves_iter = current_path[-1].add_moves
                current_sets = []

                if current_add_moves_iter != set():
                    # print (current_add_moves_iter)

                    for it in self.available:

                        # print (self.available[it])
                        # print (current_add_moves_iter in self.available[it])
                        # print (current_add_moves_iter.issubset(self.available[it]))
                        if (current_add_moves_iter.issubset(self.available[it])):
                            current_sets.append(it)
                
                if current_bp == 0 or not mode:
                    current_sets = list(range(len(destination_list)))
                    #[0,1,2,3,4]

                # if mode: 
                #     current_sets = list(range(len(destination_list)))

                # if current_bp == 0 and mode: print ("current sets:", current_sets)
                # if mode: print_d ("current sets:", current_sets)

                # "try_moves"
                for i, j in self.find_moves(current_p_table, current_sets, current_start_pt, mode):

                    # if (mode): print (i,j)

                    # if current_bp==1: print_d("i", current_string, i,j)

                    current_add_moves = current_path[-1].add_moves.copy()

                    if (i, j) in current_add_moves:
                        continue  # this optional move is already in the path

                    if (i, j) in self.moves_add:
                        current_add_moves.add((i, j))

                    # if current_bp==1:
                    #     print_d(i,j)

                    # next energy calculations
                    next_e = self.fc.eval_move(
                        current_string, i, j) + current_e
                    next_e = round(next_e, 2)


                    next_p_table = current_p_table.copy()
                    if i < 0:
                        next_p_table[-i] = 0
                        next_p_table[-j] = 0
                    else:
                        next_p_table[i] = j
                        next_p_table[j] = i

                    # next saddle energy
                    next_s = round(max(current_s, next_e), 2)

                    # if this move is valid... append to list
                    if next_s <= max_energy:
                    # if next_s < max_energy:

                        next_moves = current_moves.copy()
                        next_moves.append((i, j, next_e))

                        next_energies = current_energies.copy()
                        if next_e < s1_en:
                            next_energies.append(round(next_e-s1_en, 2))
                        else:
                            next_energies.append(0)

                        # unused?
                        en_moves = [x[2] for x in next_moves]

                        new_intermediate = Intermediate(p_table=next_p_table, mode=mode, saddle_e=next_s, current_e=next_e,
                                                        moves=next_moves, energies=next_energies, opt=[], add_moves=current_add_moves, start_pt=current_path[-1].start_pt, embedded=current_path[-1].embedded.copy())

                        new_path = current_path.copy() + [new_intermediate]
                        collect_paths.append(new_path)

            if len(collect_paths) == 0:
                break

            # first sorting step
            collect_paths.sort(key=lambda x: (x[-1].p_table, x[-1].saddle_e))

            last_ptable = []

            last_index = 0

            # for c in collect_paths:
            #     cs = RNA.db_from_ptable(c[-1].start_pt)
            # print (current_bp, cs)

            print_d(current_bp, len(collect_paths))

            # last_ptable = collect_paths[-1][0].p_table

            # last_ptable = collect_paths[0][-1].p_table
            last_start_pt = collect_paths[0][-1].start_pt

            # print_d("sort done", last_ptable, init_intermediate.p_table)

            # remove duplicates ptables
            # if current_bp+1 != current_bp_end:
            if mode:
                for i in range(len(collect_paths)):
                    c_ptable = collect_paths[i][-1].p_table
                    c_start_pt = collect_paths[i][-1].start_pt
                    if c_ptable == last_ptable:
                        # found a duplicate
                        # set saddle energy high
                        collect_paths[i][-1].saddle_e = 999
                    else:
                        last_ptable = c_ptable
                        last_start_pt = c_start_pt
            
            if not mode and current_bp+1 != self.max_bp_dist:
                for i in range(len(collect_paths)):
                    c_ptable = collect_paths[i][-1].p_table
                    c_start_pt = collect_paths[i][-1].start_pt
                    if c_ptable == last_ptable:
                        # found a duplicate
                        # set saddle energy high
                        if last_start_pt != c_start_pt and last_ptable!=[]:
                            # print ("found", i, RNA.db_from_ptable(collect_paths[i][-1].p_table), collect_paths[i][-1].moves)
                            
                            # the structure id already has an embedding
                            if c_start_pt in collect_paths[last_index][-1].embedded:
                                if collect_paths[i][-1].saddle_e <= collect_paths[last_index][-1].embedded[c_start_pt]:

                                    # print ("update emb", collect_paths[last_index][-1].embedded[c_start_pt], "to", collect_paths[i][-1].saddle_e)

                                    collect_paths[last_index][-1].embedded[c_start_pt] = collect_paths[i][-1].saddle_e
                                # else:
                                    # print ("ignore", collect_paths[i][-1].saddle_e, i)
                            else:
                                # no prev embedding (for this structure)
                                # print ("found", i, collect_paths[last_index][-1].saddle_e, last_start_pt, c_start_pt, collect_paths[i][-1].saddle_e)                                
                                # print ("start:", self.destinations[c_start_pt][0])
                                # print ("now:  ", RNA.db_from_ptable(c_ptable))
                                # print ("now:  ", RNA.db_from_ptable(last_ptable))                                
                                
                                collect_paths[last_index][-1].embedded[c_start_pt] = collect_paths[i][-1].saddle_e
                            
                            # print ("comp", last_index, collect_paths[last_index][-1].saddle_e, "i:", i, collect_paths[i][-1].saddle_e, collect_paths[last_index][-1].embedded)
                            # new: switch to best saddle en from all embedded ones


                            # for key in collect_paths[last_index][-1].embedded:
                            #     print (key, "new:", collect_paths[last_index][-1].embedded[key], "old", collect_paths[last_index][-1].saddle_e)

                            #     if collect_paths[last_index][-1].saddle_e > collect_paths[last_index][-1].embedded[key]:
                            #         print ("sw")                                    
                            #         # collect_paths[last_index][-1].embedded[last_start_pt] = collect_paths[last_index][-1].saddle_e
                            #         collect_paths[last_index][-1].start_pt = collect_paths[i][-1].start_pt
                            #         collect_paths[last_index][-1].saddle_e = collect_paths[i][-1].saddle_e
                            #         collect_paths[last_index][-1].current_e = collect_paths[i][-1].current_e




                            # continue

                            
                        
                        collect_paths[i][-1].saddle_e = 999
                    else:
                        last_index = i
                        last_ptable = c_ptable
                        last_start_pt = c_start_pt

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
            
            # debug
            # if current_bp+2 >= self.max_bp_dist:
            #     for p in collect_paths[:25]:               
            #         a = p[-1].saddle_e
            #         cs = RNA.db_from_ptable(p[-1].p_table)
            #         em = p[-1].embedded
            #         st = p[-1].start_pt
            #         print (cs, a,st, em)

            # for p in collect_paths[:5]:               
            #         a = p[-1].saddle_e
            #         cs = RNA.db_from_ptable(p[-1].p_table)
            #         em = p[-1].embedded
            #         st = p[-1].start_pt
            #         print (cs, a,st, em)



            # return valid paths if we're past the minimum bp_dist
            # if current_bp >= current_bp_end-1:  # or current_bp==30:
            if current_bp >= self.min_bp_dist:

                for p in collect_paths[:4]:               
                    a = p[-1].saddle_e
                    cs = RNA.db_from_ptable(p[-1].p_table)
                    em = p[-1].embedded
                    st = p[-1].start_pt

                    goal_dist = RNA.bp_distance(cs, destination_list[1])


                    print (cs[0:50], a,st, em, "goal d:", goal_dist)

                for i in range(len(collect_paths)):
                    end_s = collect_paths[i][-1].p_table
                    end_s = RNA.db_from_ptable(end_s)

                    if mode and end_s in destination_list:
                        yield collect_paths[i][-1]
                    if not mode and end_s in s2:
                        yield collect_paths[i][-1]

                    # sen = collect_paths[i][-1].saddle_e
                    # if sen<-15:
                    #     print ("yield", end_s)
                    #     yield collect_paths[i][-1]

                    # if collect_paths[i][-1].p_table == list(end_p_table):
                        # print("y")
                    # yield collect_paths[i][-1]

            # next iteration
            paths = collect_paths
            current_bp += 1

        # return remaining paths
        if paths:
            for path in paths:
                if path:
                    end_s = path[-1].p_table
                    end_s = RNA.db_from_ptable(end_s)

                    # print ("yield", end_s)
                    # yield path[-1]
                    if mode and end_s in destination_list:
                        yield path[-1]
                    if not mode and end_s in s2:
                        yield path[-1]


def find_path(sequence, start_list, destination_list, max_energy=float("inf"), results=1, search_width=1000, Debug=False, Verbose=False):
    """
    indirect findpath, main function

    settings:
    indirect_iterations, default value 2 (means 1 direct pass, 1 indirect pass)

    """




    if Debug:
        coloredlogs.DEFAULT_LOG_FORMAT = '%(levelname)s %(message)s'
        coloredlogs.install(level='DEBUG')


    s1 = start_list
    
    s2 = destination_list[0]

    add_moves = []

    fp_class = Fp_class(sequence, s1, s2, add_moves)


    current_width = search_width  # for the time being, only 1 pass




    # generate move list
    


    fp_class.destinations = [[i,0] for i in destination_list]

    fp_class.min_bp_dist = float("inf")
    fp_class.max_bp_dist = float("-inf")

    print_d (s1)
    print_d ("destinations:") 

    move_list = set()

    fp_class.available = dict()
    all_available = [set()] * len(destination_list)

    for a, s in enumerate(destination_list):
        bp_dist = RNA.bp_distance(s, s1)
        fp_class.destinations[a] = (s, bp_dist)

    fp_class.destinations.sort(key=lambda x: -x[1])

    for a, (s,bp_dist) in enumerate(fp_class.destinations):

        current_move_list = set()

        if bp_dist < fp_class.min_bp_dist:
            fp_class.min_bp_dist = bp_dist

        if bp_dist > fp_class.max_bp_dist:
            fp_class.max_bp_dist = bp_dist

        c_pt = fp_class.p_tables[s]
        pt1 = fp_class.p_tables[s1]

        for i in range(1, c_pt[0]):
            if pt1[i] != c_pt[i]:
                if i < pt1[i]:
                    # bp deletion
                    move_list.add((-i, -pt1[i]))
                    current_move_list.add((-i, -pt1[i]))

                if i < c_pt[i]:
                    # add bp
                    move_list.add((i, c_pt[i]))
                    current_move_list.add((i, c_pt[i]))


        fp_class.available[a] = current_move_list
        all_available[a] = current_move_list.copy()

        print_d ("dest:", a, s, bp_dist)
        print_d (len(current_move_list), current_move_list)
        # break

    destination_list = [i[0] for i in fp_class.destinations]

    #all of them combined
    fp_class.moves_add = move_list

    # One-Liner to intersect a list of sets
  
    intersection = all_available[0].intersection(*all_available)
    additional_moves = set()


    for a in all_available:
        for current_set in a:

            if current_set not in intersection:

                # add inverse move
                a = -current_set[0]
                b = -current_set[1]
                additional_moves.add((a,b)) # inverse, undo
                additional_moves.add((-a,-b))

            
    intersection = intersection.union(additional_moves)
    # print(intersection)
    # print ("add moves:", additional_moves)

    fp_class.moves_add = intersection.copy()
    # add inverse moves to undo bs...

    for k, (s,bp_dist) in enumerate(fp_class.destinations):
        for a,b in additional_moves:
            if (-a,-b) in fp_class.available[k]:                
                fp_class.available[k].add((a,b))
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

    print_d ("~FWD~~~~~~")
    result_dict = defaultdict(lambda: 9999)

    #fwd
    path_generator = fp_class.find_path_once(
        start_list, destination_list, max_energy, current_width, mode=True, sort_min=False, Debug=Debug, Verbose=Verbose)

    for path in path_generator:
        max_en = path.saddle_e

        end_s = RNA.db_from_ptable(path.p_table) #?

        # if end_s not in destination_list:
            # continue

        start_index = path.start_pt
        start_s = fp_class.destinations[start_index][0]

        print_moves(sequence, s1, start_s, path.moves)


        if max_en < result_dict[end_s]:
            result_dict[end_s] = max_en

    for key in result_dict:
        print("fwd:", key, result_dict[key]) 
# # 

    if Debug:
        coloredlogs.DEFAULT_LOG_FORMAT = '%(levelname)s %(message)s'
        coloredlogs.install(level='DEBUG')


    print_d ("~BWD~~~~~~")
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
    s1       = '.(((....))).(((((.(((..(((..(((...)))....)))...)))))))).............(((((((((((.((((...)))).))).))))))))...(((((....)))))...((((((....))))))..........'
    s2       = '.(((....))).((((.((((..(((((.((.....)))).)))...))))))))........((...(((((((((((.((((...)))).))).))))))))...(((((....)))))...((((((....))))))........))'
    intermediates = ['.(((....))).((((.((((..(((..(((...)))....)))...))))))))...........(((((((((((((.((((...)))).))).)))))))....(((((....))))))))((((((....))))))..........', '.(((....))).((((.((((..(((...((.....))...)))...))))))))...........(((((((((((((.((((...)))).))).)))))))....(((((....))))))))((((((....))))))..........', '.(((....))).((((.((((..(((.....((....))..)))...)))))))).............(((((((((((.((((...)))).))).))))))))...(((((....)))))...((((((....))))))..........']

    #15
    sequence = 'AAAACGGCUAGACGUUUAUGACUGGCGAUAAUUAUUUGCAUGAUGGAUCCAAUCCAAUGCAGUCGAGUGAAGUACCGUCGUAGUAGCUAUAACCACCUUAUUCUUGCCCCCUGCAAGCGCUAGAUACUGAUCCAUGCUAAUAAAGUAGUC'
    s1       = '.....(((((.....(((((((.((.....(((..((((((..(((((...)))))))))))...)))......))))))))))))))...((.((.((((((((((.....))))).((..(((....)))...)).))))).)).)).'
    s2       = '.(((((......)))))....((((((........((((((..(((((...)))))))))))....(((.(((..(......)..))).....)))......(((((.....))))))))))).....(((...((((.....)))))))'
    intermediates = ['.(((((......)))))((..((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).)))))))).........((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).))))))..(((.......(((.....)))))).', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).))))))....((...))((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).)))))).....(((...((((.....)))))))', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).))))))...........((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..(((((.((....((.((((...)))).))...)))))))((((.....)))).))))))...........((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((.(((....((.((((...)))).))...)))))))((((.....)))).))))))...........((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))...(((((.((...((.((((...)))).)).)).)))))(((((.....)))))))))))...........((((.....))))...', '.(((((......))))).....(((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).)))))(((....)))..((((.....))))...']
    max_energy = -12.9
    max_energy = -14.0

    # 150 nt, -20.10 kcal/mol direct, indirect -24.50
    # sequence = 'AAAAUAAUGUACCGGACAUUCGCGCACGACCACCAUAUGGCAGAGCAUGUGUCUGUGGACCCACUAUAGCUGGGGCGCUUAACCCCAGAAAAGUAUCUUCGGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAA'
    # s1       = '.............((((((..((...(..(((.....)))..).))..))))))(((((((....(((.((((((.......)))))).....)))....)))))))(((......((.((.........)).))......)))......'
    # s2       = '......(((.((((((...(((....))).(((((((((......))))))...)))((...(((....((((((.......))))))...))).)))))))).)))(((......((.((.........)).))......)))......'
    # intermediates = ['......(((.((((((.(((((....)))((((((((((......))))))...))))....(((....((((((.......))))))...))))).)))))).)))(((......((.((.........)).))......)))......', '......(((.((((((...(((....)))((((((((((......))))))...))))....(((....((((((.......))))))...)))...)))))).)))(((......((.((.........)).))......)))......', '......(((.((((((.............((((((((((......))))))...))))....(((....((((((.......))))))...)))...)))))).)))(((......((.((.........)).))......)))......', '..........((((((...(((....)))((((((((((......))))))...))))....(((....((((((.......))))))...)))...))))))....(((......((.((.........)).))......)))......']


    sequence = 'AUCAUUGCCUAAGGGACAUUUCCACGUGUAAGAAACCCGGAAACAUCGUUUGGCUUAGAAUGAACAGAUAGAUUGCCUACAACAAAAUACUGAUGACGGGGGAGCCGUGACAUGGCGGCCACCGCUAGUCGGCCUACGCUUUAAAACAUA'
    s1       = '.((((.(((.(((.((..(((((..((.......))..)))))..)).))))))..........(((.((..((..........)).)))))))))((..((.((((.....((((((...))))))..)))))).))............'
    s2       = '.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).(((.....(((.......)))....)))..((((.((..(((((......)))))..)).)..)))(((....)))..........'


    intermediates = [
        '.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
        '.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).........(((....))).......(((((..(((((..(((((...)))))..)).)))...)))))..................',
        '.((((((((.(((.((..(((((..(.((.....))).)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
        '.((((((((.(((.((..(((((...............)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
        '.((((((((.(((.((..(((((...............)))))..)).))))))....)))))..........................(((((..(((((..(((((...)))))..)).)))...)))))..................',
        '......(((.(((.((..(((((..((.......))..)))))..)).))))))..............(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................'
    ]



    # sequence = 'CACCGAGCCAGCGAGGGUCAACGCUUCACACCAUCGGGCUUGAAAGUCUUGGACCUCAGGGUAUCCGGGCCUGUAAUCACGGCUAGUGCAAUGCGGGUGUUUCGAGGGCCAGGGCUUGUUCGUUGGUCUGGUCGCCAGAGUUUGUAGAAG'
    # s2 = '...((((((..(((.(((...........))).))))))))).(((((((((.((((....((((((((((.((....)))))).........))))))....)))).))))))))).(((...(.(((((...))))).).....))).'
    # s1 = '..(((.(((((((((..((((.((((.........))))))))(((((((((.((((.((....))(.(((((((.((((.....))).).))))))).)...)))).))))))))).))))))))).)))((....))...........'
    # intermediates = ['...((((((..(((.(((...........))).))))))))).(((((((((.((((.(((((((((((((((....)).)))).........))))))))).)))).))))))))).(((...(.(((((...))))).).....))).', '...((((((..(((.(((...........))).))))))))).(((((((((.((((.(((((((((((((((....)).)))).........))))))))).)))).))))))))).......(.(((((...))))).).........', '...((((((..(((.(((...........))).))))))))).(((((((((.((((.(((((((((((((.((....)))))).........))))))))).)))).))))))))).(((...(.(((((...))))).).....))).',
    #                  '...((((((..(((.(((...........))).))))))))).(((((((((.((((.(((((((((((((.((....)))))).........))))))))).)))).))))))))).......(.(((((...))))).).........', '...((((((..(((.(((...........))).))))))))).(((((((((.((((.(((((((((((((.........)))).........))))))))).)))).))))))))).(((...(.(((((...))))).).....))).', '...((((((..(((.(((...........))).))))))))).(((((((((.((((....((((((((((((....)).)))).........))))))....)))).))))))))).......(.(((((...))))).).........']

    # find max_en between s1 and s2

    sw = 2
    max_en = pathfinder.pathfinder(
        sequence, s1, s2, search_width=RNA.bp_distance(s1, s2)*sw).max_en
    print("max_en", max_en)

    sw = 5
    search_width = int(RNA.bp_distance(s1, s2)*sw)

    # search_width = 198

    Verbose = True
    Debug = True
    # Debug = False

    # max_energy = -13.5
    # max_energy = 999

    intermediates = list(set(intermediates))

    en = find_path(sequence, s1, intermediates, max_energy=max_energy,
                   search_width=search_width, Debug=Debug, Verbose=Verbose)
    print("to s1", en)







    search_width = 500
    # search_width = 20



    print('orig findpath:')

    Verbose = False
    # Verbose = True

    destinations1 = [[i,0] for i in intermediates]
    for a, s in enumerate(intermediates):
        bp_dist = RNA.bp_distance(s, s1)
        destinations1[a] = (s, bp_dist)
    destinations1.sort(key=lambda x: -x[1])

    destinations2 = [[i,0] for i in intermediates]
    for a, s in enumerate(intermediates):
        bp_dist = RNA.bp_distance(s, s2)
        destinations2[a] = (s, bp_dist)
    destinations2.sort(key=lambda x: -x[1])

    for number, (i,_) in enumerate(destinations1.copy()):
        print ("s1: ~~~~~~~~~~")
        # print (number, i, _)
        r1 = pathfinder.pathfinder(
            sequence, s1, i, search_width=search_width, verbose=Verbose)        
        print(number, i, r1.max_en, RNA.bp_distance(s1,i))


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
