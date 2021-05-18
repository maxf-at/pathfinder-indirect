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
    add_moves:    list


def get_neighbors(fc, db=None, pt=None):
    """
    """
    if pt is None:
        pt = RNA.ptable(db)
    else:
        assert db == None
    nbrs = []
    for move in fc.neighbors(pt):        
        npt = list(pt)
        if move.is_removal():
            npt[-move.pos_3] = 0
            npt[-move.pos_5] = 0
            ndb = RNA.db_from_ptable(npt)
        elif move.is_insertion():
            npt[move.pos_3] = move.pos_5
            npt[move.pos_5] = move.pos_3
        else:
            rlog.warning(f"Are you using shift moves?")
            rlog.warning(f" shift = {move.is_shift()}")
            rlog.warning(f" pos3 = {move.pos_3}")
            rlog.warning(f" pos5 = {move.pos_5}")
            raise NotImplementedError('Are you using shift moves?')
        dG = fc.eval_move_pt(pt, move.pos_5, move.pos_3)
        if db:
            nbrs.append([RNA.db_from_ptable(npt), dG, move.pos_5, move.pos_3])
        else:
            nbrs.append([npt, dG, move.pos_5, move.pos_3])
    return nbrs


class Fp_class():
    def __init__(self, sequence, s1, s2, add_moves=[]):

        self.fc = RNA.fold_compound(sequence)
        self.sequence = sequence
        self.s1 = s1
        self.s2 = s2

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

    def find_stack(self, path, moves, max_pos):
        """
        find compensation stacks to lower the energy of e.g. saddle points
        """

        # New best result:  -5.10 kcal/mol | B:   8.50 kcal/mol | E[start]:-13.60 E[end]:-10.70 | 
        # additional moves: [(16, 24), (17, 23), (15, 25), (12, 28), (11, 29), (13, 27)]

        # common shared pairing table between s1 and s2
        p_common = list(self.p_tables[self.s1])

        for i in range(len(self.sequence)):            
            p1 = self.p_tables[self.s1][i+1]
            p2 = self.p_tables[self.s2][i+1]
            if p1 != p2:
                p_common[i+1] = 0
            # print (p1, p2)

        # print ("get nb")
        
        # print ("get nb done")    

        pt_start = [p_common]

        print ("common:")
        print (RNA.db_from_ptable(p_common))
        print ("saddle structure:", path[max_pos])

        if True:
            path[max_pos+1][0] = 80
            pt_start.append(path[max_pos+1].copy())
            path[max_pos-1][0] = 80
            pt_start.append(path[max_pos-1].copy())

        print (p_to_s(path[max_pos]))
        print ("s1s2")
        print (self.s1)
        print (self.s2)
        # print (p_common)

        # p_common = pt_start[1]

        # energy of current ptable
        base_energy = self.fc.eval_structure_pt(p_common)/100

        nbs = get_neighbors(self.fc, db=None, pt=p_common)




        

        # yield [(16, 24), (17, 23), (15, 25), (12, 28), (11, 29), (13, 27)]
        # return


        candidates = []

        # we have "moves", even without calling findpath once. 

        # this also deletes bps, yay!
        for nb, dG, i, j in nbs:
            # dont consider candidates which are in the primary move set

            if (i, j) in moves:
                continue
            if (-i, -j) in moves:
                continue

            # if current_loop[i] != current_loop[j]:
            #     continue

            current_en = round(base_energy+dG/100,3)
            # t = self.fc.eval_structure_pt(nb)/100
            candidates.append((i,j))

            # break
            # print (RNA.db_from_ptable(nb), current_en, (i, j))

        
 
 

        # allow extension of this and that:
        # either extend common ptable or one of the move_ptables in
        # the initial iteration

        current_ptable = p_common
        # current_ptable = pt_start[1]

        current_loop = self.l_tables(current_ptable)     

        def moves_to_en(moves, Verbose=False):
            # Verbose = True
            # print ("calculate moves", moves)
            
            new_ptable = current_ptable.copy()
            e = self.fc.eval_structure_pt(new_ptable)/100

            for i, j in moves:
                if i<0:
                    new_ptable[-i] = 0
                    new_ptable[-j] = 0
                else:
                    new_ptable[i] = j
                    new_ptable[j] = i

                # print(self.fc.eval_move_pt(new_ptable, i, j))

            delta_e = self.fc.eval_structure_pt(new_ptable)/100 - e
            
            # if Verbose:
            #     print(p_to_s(new_ptable), end=" ")
            #     print(delta_e, self.fc.eval_structure_pt(new_ptable)/100, moves)


            # this should use eval move step by step
            return delta_e, new_ptable


        # test = [(49,   56), (50,   55 )]
        # print("cp:", moves_to_en([]))
        # print("test:", moves_to_en(test))


        # this is stupid
        combinations = [
            [moves_to_en([(x[0], x[1])]), [(x[0], x[1])]] for x in candidates]

        all_candidates = []        


        iteration = 0
        while iteration < 10:
            iteration += 1
            # print (test)
            all_candidates += combinations
            next_combinations = []
            for current_energy, candidate_list in combinations:
                for next_candidate in candidates:
                    next_candidate_list = candidate_list.copy()

                    if next_candidate in candidate_list:
                        continue
                    start = candidate_list[-1][0]
                    end = candidate_list[-1][1]
                    # next candidate has to fit into the last one
                    if next_candidate[0] <= start or next_candidate[1] >= end:
                        continue
                    next_candidate_list += [next_candidate]
                    next_combinations.append(
                        [moves_to_en(next_candidate_list), next_candidate_list])

            if next_combinations == []:
                break
            next_combinations.sort()
            
            # keep x best candidates for next iteration
            combinations = next_combinations[0:40]

        all_candidates.sort()

        print ("bfs done")

        processed = []


        for en, candidate in all_candidates[:40]:

            

            # check if candidate is a subset of an already processed candidate
            # .((.(((...)))))(((((.....)))))...............
            # .((.((.....))))(((((.....)))))............... <-- ignore this one

            ignore = False
            for processed_candidate in processed:
                # print (processed)
                if set(candidate).issubset(processed_candidate):
                    # print ("ignore", candidate)
                    ignore = True

            if ignore:
                continue

            processed.append(set(candidate))

            _,pt_new = moves_to_en(candidate, Verbose=True)

            print(f"'{RNA.db_from_ptable(pt_new)}',")

            # print ("yield", set(candidate))

            yield candidate







        # yield [(16, 24), (17, 23), (15, 25), (12, 28), (11, 29), (13, 27)]
        # return


        # scoring system - pick the best among the additional moves
        # 






        sequence = self.sequence

        # if max_pos == 0: return set()
        # if max_pos == len(path): return set()

        for max_pos in range(len(path)):

            # print ("start find stack", max_pos, len(path), moves)

            if max_pos <= 1:
                A = path[0]
                B = path[2]
                current_ptable = path[1]
            elif max_pos >= len(path)-1:
                A = path[-3]
                B = path[-1]
                current_ptable = path[-2]
            else:
                A = path[max_pos-1]
                B = path[max_pos+1]
                current_ptable = path[max_pos]

            # workaround...
            A[0] = len(A)-1
            B[0] = len(B)-1
            current_ptable[0] = len(A)-1

            A_str = p_to_s(A)
            B_str = p_to_s(B)
            C_str = p_to_s(current_ptable)

            # quick workaround - ignore previously computed stacks
            if (A_str, B_str, C_str) in self.processed_stacks:
                continue
            self.processed_stacks.add((A_str, B_str, C_str))

            # energy of current ptable
            base_energy = self.fc.eval_structure_pt(current_ptable)/100

            # compatible = A.copy()
            # offset_loop_table = A.copy()

            loop_A = self.l_tables(A)
            loop_B = self.l_tables(B)
            current_loop = self.l_tables(current_ptable)

            # for i in range(1, len(A)):
            # if loop_A[i] == loop_B[i] and A[i]==0 and B[i]==0:
            # offset_loop_table[i] = loop_A[i] - loop_B[i]
            # if loop_A[i] == loop_B[i] and current_ptable[i] == 0:
            #     compatible[i] = True
            # else:
            #     compatible[i] = False

            allowed = ["AU", "UA", "GC", "CG", "GU", "UG"]

            # generate indirect candidates for the current structure
            candidates = []
            for i in range(1, len(A)):
                # if not compatible[i]: continue
                # offset = offset_loop_table[i]

                if current_ptable[i] != 0:
                    continue

                for j in range(i, len(A)):
                    # requirements: minimum length between i and j
                    if i+4 > j:
                        continue
                    if current_ptable[j] != 0:
                        continue

                    # if not compatible[j]: continue
                    # if offset_loop_table[j]-offset != 0: continue
                    # if current_loop[i] != current_loop[j]: continue

                    # incompatible loops (with respect to bordering structures)
                    # if loop_A[i] != loop_A[j]:
                    #     continue
                    # if loop_B[i] != loop_B[j]:
                    #     continue

                    # i and j have to be on the same loop
                    if current_loop[i] != current_loop[j]:
                        continue

                    # only consider allowed base pairs
                    seq_i = sequence[i-1]
                    seq_j = sequence[j-1]
                    if seq_i+seq_j not in allowed:
                        continue

                    # dont consider candidates which are in the primary move set
                    if (i, j) in moves:
                        continue
                    if (-i, -j) in moves:
                        continue

                    candidates.append((i, j))

            def moves_to_en(moves, Verbose=False):
                # Verbose = True
                new_ptable = current_ptable.copy()
                for i, j in moves:
                    new_ptable[i] = j
                    new_ptable[j] = i
                if Verbose:
                    print(p_to_s(new_ptable), end=" ")
                    print(self.fc.eval_structure_pt(new_ptable)/100)
                return self.fc.eval_structure_pt(new_ptable)/100

            # if max_pos == 13:
            #     print (max_pos, p_to_s(current_ptable))
            #     for candidate in candidates:
            #         moves_to_en([candidate], Verbose=True)

            # candidates_with_energies = [[moves_to_en([(x[0], x[1])]), [(x[0], x[1])]] for x in candidates]
            # candidates_with_energies.sort(key=lambda x: x[0])

            # input is a list of tuples which can potentially added [(1, 63), (1, 65), (1, 67), (1, 69), (7, 11), (7, 12)]
            # lets find out which tuple combinations have potential to lower the energy
            combinations = [
                [moves_to_en([(x[0], x[1])]), [(x[0], x[1])]] for x in candidates]
            all_candidates = []

            # combinations.sort(key=lambda x: x[0])
            # combinations = combinations[:50]

            # try: rnasubopt for candidate generation
            # cmd = f'printf "{sequence}\n{C_str}\n" | RNAsubopt -C −−enforceConstraint -e 2'
            # result = subprocess.check_output(cmd, shell=True, encoding="utf8")
            # subopt_list = result.split("\n")
            # subopt_structure = subopt_list[1].split()[0]
            # subopt_ptable = RNA.ptable_from_string(subopt_structure)
            # candidates = []
            # for i in range(1, len(current_ptable)):
            #     if current_ptable[i] == subopt_ptable[i]: continue
            #     if subopt_ptable[i] == 0: continue
            #     j = subopt_ptable[i]
            #     if i>j: continue # closing bracket

            #     if current_loop[i] != current_loop[j]: continue
            #     if (i,j) in moves: continue
            #     if (-i,-j) in moves: continue
            #     candidates.append((i, j))

            # combinations = [[moves_to_en([(x[0],x[1])]),[(x[0],x[1])]] for x in candidates]
            # all_candidates = []

            iteration = 0
            while iteration < 10:
                iteration += 1
                # print (test)
                all_candidates += combinations
                next_combinations = []
                for current_energy, candidate_list in combinations:
                    for next_candidate in candidates:
                        next_candidate_list = candidate_list.copy()

                        if next_candidate in candidate_list:
                            continue
                        start = candidate_list[-1][0]
                        end = candidate_list[-1][1]
                        # next candidate has to fit into the last one
                        if next_candidate[0] <= start or next_candidate[1] >= end:
                            continue
                        next_candidate_list += [next_candidate]
                        next_combinations.append(
                            [moves_to_en(next_candidate_list), next_candidate_list])

                if next_combinations == []:
                    break
                next_combinations.sort()
                combinations = next_combinations[0:20]

            all_candidates.sort()

            # print ("bfs done")

            # for en, candidate in all_candidates[:4]:
            # print (max_pos, candidate, en, base_energy, candidate)
            # if en<base_energy:
            #     print (max_pos, candidate, en, base_energy, candidate)

            # if max_pos == 13:
            #     # print (candidates)
            # for en, candidate in all_candidates[:10]:
                # moves_to_en(candidate, Verbose=True)
                # yield set(candidate)
                # print (max_pos, candidate, en, base_energy, candidate)

            if all_candidates != [] and all_candidates[0][0] < base_energy:
                
                # print (p_to_s(current_ptable), "+", all_candidates[0][-1])
                moves_to_en(all_candidates[0][-1], Verbose=True)

                # for move in 


                # print(f"'{RNA.db_from_ptable(pt_new)}',")
                # print ("c:", all_candidates[1][-1], moves_to_en(all_candidates[1][-1], Verbose=True))
                # print ()
                # i = all_candidates[0][1][0][0]
                # j = all_candidates[0][1][0][1]
                # print ("candidate:", all_candidates[0][0]-base_energy, all_candidates[0][1], (i,j) in moves, (-i,-j) in moves)
                yield set(all_candidates[0][-1])

            # if all_candidates != [] and len(all_candidates)>1 and all_candidates[1][0]<base_energy:
                # # print ("c:", all_candidates[1][-1], moves_to_en(all_candidates[1][-1], Verbose=True))
                # yield set(all_candidates[1][-1])

        return

    def find_moves(self, s_ptable, t_ptable):
        """
        generator function, yields possible structures 1 move away
        from the original structure by finding fitting i and j with
        RNA pair and loop tables
        s_ptable: current ptable
        t_ptable: s2 end ptable
        """

        # loop table
        ls = self.l_tables(s_ptable)

        # extra indirect moves: i and j have to be unpaired, and on the same loop
        for i, j in self.moves_add:
            if i > 0 and s_ptable[i] == 0 and s_ptable[j] == 0 and ls[i] == ls[j]:
                # self.moves_add.remove((i,j))
                # print ("try", i, j)
                yield i, j
                # pass
            if i<0:
            # if i<0 and -j == s_ptable[-i]:
                # and -i==s_ptable[-i]:
                # print ("try -", i, j, s_ptable[i], s_ptable[-i])
                # self.moves_add.remove((i,j))
                yield i,j

        for i in range(len(s_ptable)):
            if i == 0:
                continue

            if s_ptable[i] == 0 and t_ptable[i] > i:
                j = t_ptable[i]
                # found j has to be empty and currently on the same loop as i
                if s_ptable[j] == 0 and ls[i] == ls[j]:
                    yield i, j
            # test for bp removal: i has to be paired with a different j in s2
            j = s_ptable[i]
            # dont remove things which are present in s2
            if s_ptable[i] > i and s_ptable[i] != s_ptable[j] and\
                    s_ptable[i] != t_ptable[i] and s_ptable[j] != t_ptable[j]:
                yield -i, -j

    def find_path_once(self, s1, s2, max_energy, width, mode=True, sort_min=False, Debug=False, Verbose=False):
        """
        main findpath algorithm (bounded BFS)
        """

        if Debug:
            coloredlogs.DEFAULT_LOG_FORMAT = '%(levelname)s %(message)s'
            coloredlogs.install(level='DEBUG')

        # caching
        evals = self.evals
        p_tables = self.p_tables
        l_tables = self.l_tables
        bp_dist = self.bp_dist
        # s1 = fp_class.s1
        # s2 = fp_class.s2

        e1 = self.evals[s1]
        e2 = self.evals[s2]

        runtimes = 0
        current_bp = 0
        current_bp_end = bp_dist(s1, s2)      # current_bp_end = 4
        # paths = [(list(p_tables[s1]), float("-inf"), evals[s1], [])]

        s1_en = round(evals[s1], 2)
        end_p_table = self.p_tables[s2]
        init_intermediate = Intermediate(p_table=list(p_tables[s1]), mode=mode, saddle_e=float(
            "-inf"), current_e=evals[s1], moves=[(0, 0, s1_en)], energies=[0], opt=0, add_moves=[])
        # initial path start with 1 intermediate
        init_path = [init_intermediate]
        # paths start with 1 initial path
        paths = [init_path]

        # dont stop at current_bp_end, consider potential indirect moves
        while (current_bp != current_bp_end+2*len(self.moves_add)):

            # collect all new paths here (next iteration)
            collect_paths = []

            for current_path in paths:
                current_p_table = current_path[-1].p_table
                current_e = current_path[-1].current_e
                current_s = current_path[-1].saddle_e
                current_string = p_to_s(current_p_table)
                current_moves = current_path[-1].moves
                current_energies = current_path[-1].energies

                # "try_moves"
                for i, j in self.find_moves(current_p_table, end_p_table):

                    current_add_moves = current_path[-1].add_moves.copy()

                    if (i, j) in current_add_moves:
                        continue  # this optional move is already in the path

                    if (i, j) in self.moves_add:
                        current_add_moves.append((i, j))

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
                                                        moves=next_moves, energies=next_energies, opt=[], add_moves=current_add_moves)

                        new_path = current_path.copy() + [new_intermediate]
                        collect_paths.append(new_path)

            # first sorting step
            collect_paths.sort(key=lambda x: (x[-1].p_table, x[-1].saddle_e))

            last_ptable = []
            # last_ptable = collect_paths[-1][0].p_table
            # last_ptable = collect_paths[0][-1].p_table
            print_d("sort done", last_ptable, init_intermediate.p_table)

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

            # return valid paths if we're past the minimum bp_dist
            if current_bp >= current_bp_end-1:# or current_bp==30:
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


def find_path(sequence, s1, s2, add_moves=[], results=1, indirect_iterations=2, search_width=1000, Debug=False, Verbose=False):
    """
    indirect findpath, main function

    settings:
    indirect_iterations, default value 2 (means 1 direct pass, 1 indirect pass)

    """

    fp_class = Fp_class(sequence, s1, s2, add_moves)
    max_energy = float("inf")

    current_width = search_width  # for the time being, only 1 pass

    indirect_moves = [add_moves]
    best_max_en = float("inf")
    best_max_pos = 0
    best_path = False
    best_indirect_moves = []
    best_ptables = []
    best_indirect_move_count = 0

    return_value = 9999

    init_max_en = False

    for iteration in range(indirect_iterations):

        # first pass is direct, then iterate over different indirect moves

        next_indirect_moves = []

        # iterate over moves from the last iteration - this is a single iteration in the first run.
        for indirect_move in indirect_moves:

            if indirect_move == None:
                indirect_move = []  # first iteration only

            indirect_move = set(best_indirect_moves) | set(indirect_move)
            if Verbose:
                print("Iteration", iteration,
                  "launching findpath with addtional moves:", indirect_move)

            # add a set of optional indirect moves to the direct move set, e.g. add {(6, 12), (5, 13)}
            fp_class.moves_add = indirect_move

            # main findpath_once function call with supplied settings
            path_generator = fp_class.find_path_once(
                s1, s2, max_energy, current_width, mode=True, sort_min=False, Debug=Debug, Verbose=Verbose)

            # iterate over all currently found paths. somepaths might not be optimal
            for path in path_generator:
                # path is an intermediate struct, containing saddle_e, moves, etc.

                if path.saddle_e < max_energy:
                    max_energy = path.saddle_e

                e_0 = fp_class.evals[s1]
                e_1 = fp_class.evals[s2]
                current_ptable = list(fp_class.p_tables[s1])

                ptables = []  # list of ptables according to moves, should not be needed later

                # workaround convert moves to list of ptables, del this later
                current_used_indirect_moves = []  # which indirect moves are actually used
                current_moves = []
 
                saddle_pos = None

                for pos, (i, j, e) in enumerate(path.moves):

                    if e == max_energy:
                        # print ("pos", pos)                        
                        saddle_pos = pos

                    current_moves.append((i, j))
                    if (i, j) in indirect_move:
                        current_used_indirect_moves.append((i, j))
                    if i < 0:
                        current_ptable[-i] = 0
                        current_ptable[-j] = 0
                    else:
                        current_ptable[i] = j
                        current_ptable[j] = i
                    ptables.append(current_ptable.copy())



                # save best path for print output / next iteration
                if max_energy < best_max_en or \
                   (max_energy == best_max_en and len(current_used_indirect_moves) < best_indirect_move_count):
                   # or save best path with same energy, but fewer indirect moves: take this one instead

                    if not init_max_en:
                        init_max_en = max_energy

                    best_moves = current_moves.copy()
                    best_max_pos = saddle_pos

                    best_indirect_move_count = len(current_used_indirect_moves)
                    best_max_en = max_energy
                    best_path = path.moves.copy()
                    best_ptables = ptables.copy()
                    best_indirect_moves = current_used_indirect_moves.copy()
                    barrier = max_energy-e_0
                    if Verbose:
                        print(
                        f"New best result: {max_energy:6.2f} kcal/mol | B: {barrier:6.2f} kcal/mol | E[start]:{e_0:6.2f} E[end]:{e_1:6.2f} | additional moves: {current_used_indirect_moves}")

                # elif max_energy == best_max_en:
                #     print(
                #         f"Result:          {max_energy:6.2f} kcal/mol | B: {barrier:6.2f} kcal/mol | E[start]:{e_0:6.2f} E[end]:{e_1:6.2f} | additional moves: {current_used_indirect_moves}")


                # indentation below...
                # call find stack for every path? or every iteration

                # find_stack call: here we find indirect moves, which 
                # are used in the next iteration
                # if iteration+1 != indirect_iterations:  # dont find new indirect moves in the last iteration                
                #     print ("call find_stack")
                #     for current_indirect_moves in fp_class.find_stack(ptables, current_moves, best_max_pos):
                #         # print ("find move", current_indirect_moves, next_indirect_moves)
                #         if current_indirect_moves not in next_indirect_moves:
                #             next_indirect_moves.append(current_indirect_moves)

            if iteration+1 != indirect_iterations:  # dont find new indirect moves in the last iteration                
                print ("call find_stack for best path")
                for current_indirect_moves in fp_class.find_stack(best_ptables, best_moves, best_max_pos):
                    # print ("find move", current_indirect_moves)
                    if current_indirect_moves not in next_indirect_moves:
                        print ("add", current_indirect_moves)
                        next_indirect_moves.append(current_indirect_moves)





                # break # stop at 1 path

        # print path during last iteration
        if iteration+1 == indirect_iterations:
            # print(best_path)
            return_value = print_moves(sequence, s1, s2, best_path, Verbose=Verbose)
        # prepare for next indirect iteration
        else:
            # print ("found indirect moves:")
            # print (next_indirect_moves)
            indirect_moves = next_indirect_moves

        # best_path = [(0, 0), (-31, -45), (-35, -41), (-34, -42), (-33, -43),
        #  (-32, -44), (-5, -11), (-4, -12), (6, 12), (7, 11),
        #  (-3, -13), (5, 13), (-2, -14), (-1, -15), (4, 42), (3, 43), (2, 44), (1, 45),
        #  (-5, -13),(5, 41), (-6, -12), (6, 40), (-7, -11), (7, 39),
        #  (8, 38), (9, 37), (10, 36)]

    # code here to not just return nothing
    return return_value


if __name__ == '__main__':

    # various random sequences

    # # 60.1 inner
    # sequence = "GCCAACAAACCGUGAUGGGCUAUGUUUGAUUCAUCCUAUUUAUGUUUUUCGAGAUGCGCG"
    # s1       = "...........(((..(((..((.....))...)))....)))................."
    # s2       = "...........((((((((..(((.......)))))))).)))................."

    # sequence = "CACGCUGGGGAAUUCCGACGAUAAAUAACCCUUCGCGCCCAUAGGAACGUACACAAAUCUGAAGCUUACG"
    # s1       = "..(((.(((...................)))...)))((....))..((((...............))))"
    # s2       = "..(((.((((..................))))..)))................(...............)"

    # sequence = "GCCGCCUGAGCCAUUUACAAACUAUCGGAUACAGUGUGUCUUAGAUGUUGUGGCAACAUUAGAUAAGGUA"
    # s1       = "(((..(((((((((..(((..(((..((((((...))))))))).))).)))))....))))....)))."
    # s2       = "...((((..(((((..(((..(((...(((.......))).))).))).)))))...........))))."


    # sequence   = "UAGGGGUGCACUAAAGCUGGUAUCCCCUAUGAGUGGAUAAAUGAUACAGGUCACCCUACGACAUAUACGC"
    # s1         = "(((((.((.(((........((((((......).))))).........))))))))))............"
    # s2         = "..((((((..((........(((((.........)))))........))..)))))).((.......))."

    # additive example
    # sequence   = "CCAGCGUAUUAGUUAUGGCCUGGAGGUAGAAGCGUUAGAGCAAUACUUCUACAGAGACCACGUGAGGUAG"
    # s1         = "((((..((.......))..))))..(((((((.............))))))).................."
    # s2         = "((((((((.....))).)).)))..(((((((.((.......)).)))))))....(((......))).."

    # sequence   = "UCUUGAACCCAUGGCGUUUUCAACCGACAUCCUGCUCCCGCAAUCACCUAGGUUAAGGGUUCUUCAAGGA"
    # s1         = "((((((((((.((((.........(.......(((....)))........))))).))))...))))))."
    # s2         = "....((((((.(((((((...)))........(((....))).........)))).))))))........"

    # simple additive
    # sequence = "UCACUGAGGCUUGUUCGCAAAUCACUGCAAUUAGAUAUGACUCACGAUAUGGGGCACGGUGCAUACAUAC"
    # s1       = ".(((((.(.((((((((....((...............))....))).))))).).)))))........."
    # s2       = ".(((((..((((..(((..(.(((.((........))))).)..)))....)))).)))))........."

    # ~~~

    # Tabu paper toy example
    # sequence = "CGCGACGGCUACGCGACGGCAAUGCCGUUGCGAAGCCGUCGCGAUC"
    # s1 = "(((((((((..............))))))))).............."
    # s2 = "...........(((((((((..............)))))))))..."

    # rna2dfold example
    sequence = "GGGCGCGGUUCGCCCUCCGCUAAAUGCGGAAGAUAAAUUGUGUCU"
    s1 = "(((((.....)))))(((((.....)))))(((((.....)))))"
    s2 = "((((((((((.....(((((.....))))).....))))))))))"

    # EApath supplement

    # rb2: indirect move inbetween s1 / s2
    # sequence = "GGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCGUAAACUCUUGAUUAUGAAGUCUGUCGCUUUAUCCGAAAUUUUAUAAAGAGAAGACUCAUGAAU"
    # s1       = "............((((((.........))))))........((((((.......)))))).((((((.((((.((.((((((..........)))))).)))))))))))).."
    # s2       = "(((((((((...((((((.........))))))........((((((.......))))))..))))))))).........................................."

    # rb3
    # sequence = "CUCUUAUCAAGAGAGGUGGAGGGACUGGCCCGAUGAAACCCGGCAACCAGCCUUAGGGCAUGGUGCCAAUUCCUGCAGCGGUUUCGCGUUGAAAGAUGAGAGAUUCUUGUAGUCUCUUCUUUUAGCGAAGGGACUUUUUUU"
    # s1       = "((((((((....(.(((...(((.....)))......))))(((..((((((....))).))).)))........(((((......)))))...)))))))).........((((((((.......))))))))......."
    # s2       = "............(.(((...(((.....)))......))))(((..((((((....))).))).)))........(((((......)))))(((((.((((...........)))).)))))..................."

    # lots of nbs
    # rb5 -26.40 kcal/mol vs -23.60 kcal/mol direct
    # sequence = "GGGAAUAUAAUAGGAACACUCAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGACUAUGGGUGAGCAAUGGAACCGCACGUGUACGGUUUUUUGUGAUAUCAGCAUUGCUUGCUCUUUAUUUGAGCGGGCAAUGCUUUUUUUAUUCUCAUAACGGAGGUAGACAGGAUGGAUCCACUGA"
    # s1       = "................((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))........(((((........)))))............((((((((((((((.......))))))))))))))..........................................."
    # s2       = ".....................(((...(.(((((.......))))).)........((((((.......))))))..)))((((((((((((.(((((........)))))..............))))))))))))................................................................"

    # dsrA
    # sequence = "ACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUAAGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUU"
    # s1 = "..(((((((.....)))))))...(((((((((((((.......))))))).)))))).((((((((((.....))))))))))."
    # s2 = "..(((((((.....)))))))................(...(((((....)))))...)((((((((((.....))))))))))."
    # s1, s2 = s2, s1


    # HIV leader
    # sequence = "GGUCUCUCUGGUUAGACCAGAUUUGAGCCUGGGAGCUCUCUGGCUAACUGGGAACCCACUGCUUAAGCCUCAAUAAAGCUUGCCUUGAGUGCUUCAAGUAGUGUGUGCCCGUCUGUUGUGUGACUCUGGUAACUAGAGAUCCCUCAGACCCUUUUAGUCAGUGUGGAAAAUCUCUAGCAGUGGCGCCCGAACAGGGACUUGAAAGCGAAAGGGAAACCAGAGGAGCUCUCUCGACGCAGGACUCGGCUUGCUGAAGCGCGCACGGCAAGAGGCGAGGGGA"
    # s1       = "((..((.(..(((((.(((((...((((......))))))))))))))..)))..)).((.(((..(((((........(((((.((.((((((((.((((.....(((.((((((((.((..(((((((.....(((....)))...(((((((.((((.((((((.....)))).)).)))).(((.....)))..........)))))))..)))))))..)).....))))...))))..))))))))))))))).)).))))))))))))).))."
    # s2       = "((..((.(..(((((.(((((...((((......))))))))))))))..)))..)).((.(((..(((((........(((((.((.(((((((.(((((.....(.(((...(((.((((.(((((((..((.(((....))))).(((((((.((((.((((((.....)))).)).))))((((.....))).)........)))))))..)))))))(((....)))..)))).))).)))))))))))))))).)).))))))))))))).))."

    # barriers
    # ((((((((((...(.(((((.....))))).)...)))))))))) (-17.70) L0002
    # ((((((((((.....(((((.....))))).....)))))))))) (-17.70) I
    # (((((((((......(((((.....)))))......))))))))) (-16.60) I
    # ((((((((.......(((((.....))))).......)))))))) (-15.00) I
    # (((((((........(((((.....)))))........))))))) (-14.40) I
    # (((((((...)....(((((.....))))).........)))))) ( -8.90) I
    # (((((((...))...(((((.....)))))..........))))) ( -7.00) S
    # (((((((...)))..(((((.....)))))...........)))) (-10.40) L0020
    # .((((((...)))..(((((.....)))))...........))). ( -9.80) I
    # ((.((((...)))..(((((.....)))))...........))). ( -6.20) S
    # ((..(((...)))..(((((.....)))))............)). ( -6.30) L0120
    # (...(((...)))..(((((.....))))).............). ( -6.00) S
    # (...(((...))).)(((((.....)))))............... ( -7.70) I
    # ((..(((...)))))(((((.....)))))............... (-10.50) L0019
    # (((..((...)))))(((((.....)))))............... (-10.40) S
    # ((((..(...)))))(((((.....)))))............... (-11.40) I
    # (((((.....)))))(((((.....)))))............... (-16.80) L0003
    # (((((.....)))))(((((.....))))).(...........). (-12.70) S
    # (((((.....)))))(((((.....))))).((.........)). (-13.90) I
    # (((((.....)))))(((((.....))))).(((.......))). (-15.90) I
    # (((((.....)))))(((((.....))))).((((.....)))). (-17.20) I
    # (((((.....)))))(((((.....)))))(((((.....))))) (-18.10) L0001


    # sequence = "CACAACCUUACCGCUGUUCUGCACACUGGAUGAUCUCUCUGGUAGUUUACCUAACCGGAUUUUAGCUCAAGGUUAAGACG"
    # s1 = "..((.((.....((......)).....)).)).....((((((((.....))).)))))((((((((...)))))))).."
    # s2 = "...((((((...((((.(((.((...)).).))....((((((..........))))))...))))..))))))......"

    # additive example from 500 dataset
    # sequence = "GGUUGUCACCAGGACCCUAGCUCCAUGACGAGGAAGAACGGCUGCGACUUUUGGGGGAGUCCCGAGAACUCACGCUUUAA"
    # s1 = "..((((((...(((.......))).))))))........(((.(.((.(((((((.....)))))))..)).))))...."
    # s2 = "((((........))))((..(((......)))..))...(((.(.((.((((.(((....))).)))).)).))))...."

    # additive ex 2
    # sequence = "ACAAACUUUAUGAGUGGACGCAUAGUUCUGGCUUAUUCGGGGCGUUUAGCAAGGGUGACAACGCAGUUCGGUUAAUGGUC"
    # s1 = ".....((..((((((.(((.....).))..))))))..))(((..(((((..(((((......)).))).)))))..)))"
    # s2 = ".........((((((((((.....))))..))))))((((((((((.............)))))..)))))........."

    # has to be inversed
    # sequence = "CGUUGGCUGCGUCAAAUCUUAUCAGUCCUUACGACACAAAGCCUGCUGUCCGCCGUUGUUGGAGAGGCGGCCUAGAGCAG"
    # s1 = "(((.(((((.((........))))).))..))).........(((((.((((((.((.....)).)))))...).)))))"
    # s2 = "....(((((.(((...................))).)..))))((((.(..((((((........))))))..).))))."
    # add_moves = [(6,26),(5,27)]

    # 1
    # sequence = "AAGAAGACCUCAAUCGAAUCACGGGCAAGUCCGACGAGGAACGCCUAGGCGAGGUGAUCGGCCCGAUCUUAAUGUAGGAU"
    # s1 = "..((.((......))...)).(((((..((((..(((((....)))...))..).)))..)))))((((((...))))))"
    # s2 = ".......(((...........(((((...(((.....))).(((((.....)))))....)))))..........))).."
    # # add_moves = {(32, 60), (37, 55), (43, 50), (42, 51), (33, 59), (34, 58), (44, 49)}
    # # add_moves = {(43, 50), (42, 51), (44, 49)}
    # add_moves = [(35, 49), (32, 53), (51, 59), (52, 58), (50, 60), (44, 49)]

    # 2
    # sequence = "GGAAGCCGGCGAGGCAGUACCAUUAUAUAGUUUGUCUUCCAAGAAUGGGUACGACCGCGGGACCGUUCGGUUAUCGUCUG"
    # s1 = ".((((((((.(((((((..............)))))))))..((((((...((....))...)))))))))).))....."
    # s2 = "((.((((((((.((..((((((((...................))))).)))..)).......)).)))))).))....."

    # add_moves = [(40, 47), (39, 48), (41, 46), (38, 49), (13, 63)]
    # add_moves = [(40, 47), (39, 48), (41, 46), (38, 49)]
    # add_moves = [(48, 64), (49, 63), (50, 62), (30, 36), (29, 37)]

    # at higher search width, this wins (even 20)
    # add_moves = [(40, 47), (39, 48), (41, 46), (38, 49), (14, 58), (15, 57), (17, 55), (18, 54), (20, 49)]
    
    # add_moves = [(40, 47), (39, 48), (41, 46), (38, 49), (14, 58), (15, 57)]
    # add_moves = [(52, 60), (48, 56), (49, 55), (50, 54), (40, 47), (41, 46), (30, 36), (29, 37)]


    # 2 alt 1-4
    # sequence = "GGAAGCCGGCGAGGCAGUACCAUUAUAUAGUUUGUCUUCCAAGAAUGGGUACGACCGCGGGACCGUUCGGUUAUCGUCUG"
    # s1 = "((.((((((.(((((((.((.........)))))))))))..((((((...((....))...)))))))))).))....."
    # s2 = "((.((((((((.((..((((((((....((.....))......))))).)))..)).)).......)))))).))....."
    # s1, s2 = s2, s1

    # 3
    # sequence = "UGUCGGGAUAGUCGAGGGAUGAUUUCUCUUUAGUGCCGGCACUCGAGGCCUCGUAACCCGACCCAAUGGGAUCCGAAUAU"
    # s1 = "((((((.((....((.(((....)))))....)).)))))).(((.((.(((((...........))))).)))))...."
    # s2 = ".((((((......((((((...))))))...((.(((.........))))).....)))))).................."


    # 18
    # sequence = "GUGCCAUUGGGAAUCGCCCAGCGGGCGUUCAGAAUAGAAUCCCUACACUUCUAAGAGGAACUCAUGUGAUGUUUGGCUAC"
    # s1 = "..(((((((((.....))))).(((.((((......)))))))....(((....))).(((((....)).)))))))..."
    # s2 = "((((((....(((.(((((...))))))))....(((.....)))(((.....((.....))...))).....))).)))"
    # add_moves = [(51, 58), (52, 57), (50, 59), (49, 60), (10, 17), (9, 18)]

    # 27
    # sequence = "CACGACGUUAUAGCAUUCACUAAGUAUACUAUUGAUGUAGAAAUGCUAAUGCGCUGGUCAGGUGCAAUACCGGCAGCAGA"
    # s1 = "..........(((((((...((....))((((....)))).))))))).(((((((((..........)))))).))).."
    # s2 = "...(.(((..((((((((......(((((....).))))))).)))))))))((((.((.((((...))))))))))..."
    # add_moves = [(20, 40), (21, 39), (22, 38), (25, 35), (24, 36)]
    # add_moves = [(20, 37), (19, 38), (25, 31), (24, 32), (23, 33), (26, 30), (17, 40)]
    # add_moves = [(20, 37), (19, 38), (25, 31), (24, 32), (23, 33), (26, 30)]


    # 33 this ex is a good benchmark to beat...
    # sequence = "UAUUUUGUGUCACGCAAUUGGUUCACCCACUAUGCACGAAGUAUGCGAAGUACAGUCUAGUUAACAGUAAUAAUGCCUGG"
    # s1 = "((((((((((.(((((..(((.....)))...))).....))))))))))))(((.(.................).)))."
    # s2 = "((((((((.....(((..((((......))))))))))))))).((........))........(((.........)))."
    # add_moves = [(55, 65), (56, 64), (53, 67), (52, 68), (51, 69), (46, 75), (45, 76), (47, 74), (11, 38), (10, 39)]
    # add_moves = [(11,38), (10,39)]
    # add_moves = [(44, 50), (43, 51), (42, 52), (41, 53), (11, 38), (10, 39)]

    # 48
    # sequence = "GAUCAGUUGACCAAUUCCAGCACGUGAGUCUUAACGGAGUUUGCAUUCACAGUUUUAACUGGUACAGGUCGUAAAUCUGU"
    # s1 = ".(((((((((.............(((((((((....))).....))))))....)))))))))((((((.....))))))"
    # s2 = "(((....(((((....((((...(((((((............).))))))........))))....)))))...)))..."
    # add_moves = [(23, 52), (22, 53), (33, 41), (34, 40), (35, 39), (9, 18)]
    # s1, s2 = s2, s1

    # 71 hard
    # sequence = "ACUUCUCAGGGUUUGCCGGCGCCUCCCACAAAGCACGUCUCUGCGGGGUGUUGAUUAGUUCCGUCACUAUCAAAACUGCC"
    # s1 = "........(((.(((.(((((..((((.....(((......))))))))))))..))).))).................."
    # s2 = "........((.....))((((..((((.((.((.....)).)).))))..((((.((((......))))))))...))))"
    # add_moves = [(-24, -48), (22, 48), (-25, -47), (23, 47), (-21, -49), (-20, -50), (-19, -51), (-18, -52), (-17, -53), (-17, -53), (-18, -52), (-19, -51), (-20, -50), (-21, -49)]
    # s1 = "........(((.(((.(((((((..((.....((........)))))))))))..))).))).................."
    # add_moves = [(22,   50), (21,51)]
    # add_moves = [(-24, -48), (22,   48), ( 23,   47), (-25, -47), (48, 76), (49, 75)]
    

    # 75
    # sequence = "CGAUCGCCGAACUUACGUCGGGUCGGAAGAGGGUACAUGCGGGUCAAAGCUAGUGACGUGCCUACUCCCGCAUCGGUUCU"
    # s1 = "((((..(((.((....)))))))))...........(((((((.....((.........)).....)))))))......."
    # s2 = "((((((.((......)).).)))))((.(.(((.....(((.((((.......)))).))).....))).).))......"
    # add_moves = [(10, 18), (36, 75), (32, 79), (33, 78), (31, 80)]
    # s1, s2 = s2, s1


    # sequence = 'UAAAAUGAUCACGGUUUCAGCUUUGGACAGGGCGUUCCACUAAACUCCUGGUGACAUAGAUAUAUUGGAUUGCAACUACUCGUCGGUCCGGUUGGCGUUC'
    # s1       = '....(((.((((...............(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........'
    # s2       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'

    # s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'
    # s2       = '....(((.(((((((....))).....(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........'
    # add_moves = [(40, 50), (37, 53), (36, 54), (64, 77), (63, 78), (50, 55), (49, 56), (59, 76), (58, 77), (57, 78), (38, 50), (39, 49), (37, 51), (35, 53), (67, 76), (60, 66), (53, 73), (52, 74), (55, 72), (56, 71), (51, 76), (50, 77), (49, 78), (52, 75), (32, 47), (31, 48), (30, 49), (53, 73), (52, 74), (60, 66), (50, 76), (50, 55), (67, 73), (66, 74), (65, 75), (59, 81), (60, 80), (53, 84), (54, 83), (55, 82), (52, 85), (59, 76), (58, 77), (57, 78), (59, 81), (60, 80), (76, 82), (75, 83), (73, 85), (73, 82), (25, 55), (24, 56), (21, 59), (22, 58), (15, 94), (18, 25), (19, 24), (17, 26), (16, 27), (76, 82), (75, 83)]




    # sequence = 'AAGAAGACCUCAAUCGAAUCACGGGCAAGUCCGACGAGGAACGCCUAGGCGAGGUGAUCGGCCCGAUCUUAAUGUAGGAU'
    # s1       = '.....((......))......(((((..((((..(((((....)))...))..).)))..)))))((((((...))))))'
    # s2       = '.......(((...........(((((...(((.....))).(((((.....)))))....)))))..........)))..'
    # add_moves = [(35, 54), (42, 48), (34, 55), (39, 50), (42, 51), (43, 50), (44, 49), (32, 54)]

    # sequence = 'CUUCCAGAGUUUCGUUGAACCGGCAAACGGCCACUAAGGAUAGCGCCUCCCGCUUGUCCCUGGACCGGCGUGCAUCGUAU'
    # s1       = '......................(((..((((((....((((((((.....)))).)))).))).)))...))).......'
    # s2       = '..(((((.((((....)))).(((.....))).....((((((((.....))).))))))))))...(((.....)))..'
    # add_moves = [(32, 63), (34, 63), (35, 62), (36, 61), (28, 68), (30, 66), (29, 69), (31, 67), (28, 70), (27, 71), (32, 67), (31, 68), (23, 73), (24, 72), (25, 71), (30, 69), (28, 70), (27, 71), (22, 28), (21, 29), (20, 30), (7, 76), (8, 75), (13, 70), (14, 69), (9, 73), (10, 72), (15, 68), (16, 67), (17, 66)]

    # sequence = 'CUGGGUCGUGGUGGCCUUUUAGAUACGAUUCACGAACGUAGCACGUUUCGGUCUCCGGAGACGCAAUGAUCUCGAGGGUA'
    # s1       = '.(((((((((.............)))))))))(((((((...))))).)).((((.(((..........))).))))...'
    # s2       = '.(((((((((.............)))))))))(((.(((.((...((((((...))))))..)).)))...)))......'
    # add_moves = [(57, 73), (58, 72), (60, 71), (61, 70), (59, 72), (58, 73), (55, 76), (54, 77), (53, 78), (34, 48), (48, 58), (47, 59), (46, 60), (45, 62), (44, 63), (46, 61), (47, 60), (48, 59), (40, 67), (39, 68), (48, 58), (47, 59), (46, 60)]

    # sequence = 'GGAAGCCGGCGAGGCAGUACCAUUAUAUAGUUUGUCUUCCAAGAAUGGGUACGACCGCGGGACCGUUCGGUUAUCGUCUG'
    # s1       = '((.((((((((.((..((((((((....((.....))......))))).)))..)).)).......)))))).)).....'
    # s2       = '((.((((((.(((((((.((.........)))))))))))..((((((...((....))...)))))))))).)).....'

    # s2       = ".........((.((..((((((((....((.....))......))))).)))..)).))....................."
    # s2       = ".............(..((((((((....((.....))......))))).)))..).(((.((((....))))..)))..."

    # s1       = ".............(..((((((((....((.....))......))))).)))..).(((.((((....))))..)))..."
    # s2       = '((.((((((.(((((((.((.........)))))))))))..((((((...((....))...)))))))))).)).....'

    # s1       = '((..(((.....))).((((((((....((.....))......)))).))))..))((((((((....)))).))))...'                
    # s1       = '............((..((((((((......((((.....))))))))).)))..))((((((((....)))).))))...'
    # s2       = '((.((((((.(((((((.((.........)))))))))))..((((((...((....))...)))))))))).)).....'    

    # s1       = '((.((((((((.((..((((((((....((.....))......))))).)))..)).)).......)))))).)).....'
    # s2       = '...((((((.(((((((.((.((...)).)))))))))))..(((..(((.((....))..))).)))))))........'


    # add_moves = [(8, 64), (9, 63), (7, 65), (2, 68), (3, 67), (4, 66), (6, 61), (7, 60), (1, 64), (2, 63), (9, 58), (10, 57), (63, 70), (64, 69), (62, 71), (10, 59), (11, 58), (61, 72), (58, 76), (59, 75), (57, 77), (60, 74), (6, 14), (7, 13), (5, 15), (47, 56), (48, 55), (1, 21), (2, 20), (29, 37), (30, 36), (1, 40), (2, 39), (3, 38), (21, 34), (22, 33), (4, 37), (17, 28), (18, 27), (19, 26), (5, 36), (8, 40), (9, 39), (46, 57), (58, 65), (59, 64), (6, 70), (7, 69), (5, 71), (4, 72), (1, 75), (2, 74), (8, 68), (9, 67), (39, 48), (40, 47), (41, 46), (9, 66), (10, 65), (48, 64), (49, 63), (50, 62), (8, 40), (9, 39)]
    
    # -4.00
    # add_moves = [(49,63),(50,62),(48,64)]

    # works with 1200
    # add_moves = [(57, 77),(58, 76), (59, 75),  (60, 74),  (61, 72),  (62, 71), (63, 70), (64, 69),]
    # add_moves += [(47, 56), (48, 55)]
    # add_moves += [(5,15), (6,14), (7,13)]
    # add_moves += [(1,40), (2, 39), (3, 38)]
    # add_moves += [(-1, -75), (-2, -74), (-4, -72), (-5, -71), (-6, -70), (-7, -69), (-8, -68),]
    # add_moves += [( 8,  64), (9, 63)]

    # todo: Helices: group them together...


    # add_moves = []
    # add_moves = [(-9, -67), (9, 67)]
    
    
    # add_moves += [(-1, -75), (-2, -74), (-4, -72), (-5, -71), (-6, -70), (-7, -69)]
    

    # s1 = "(((((..............................)))))......((......))((((((((....)))).))))..."
    # s2       = '((.((((((.(((((((.((.........)))))))))))..((((((...((....))...)))))))))).)).....'

    
    # add_moves = [(-9, -67), (-8, -68), (8, 64), (9, 63), (-1, -75), (-2, -74), (-4, -72), (-5, -71), (-6, -70), (-7, -69), (7, 65), (2, 68), (3, 67), (4, 66), (-4, -66), (-2, -68), (-3, -67), (-7, -65), (-8, -64), (-9, -63), (6, 61), (7, 60), (1, 64), (2, 63), (-10, -59), (-11, -58), (9, 58), (10, 57), (-1, -64), (-2, -63), (63, 70), (64, 69), (62, 71), (-9, -58), (-10, -57), (10, 59), (11, 58), (-6, -61), (-7, -60), (61, 72), (-10, -59), (-11, -58), (58, 76), (59, 75), (57, 77), (60, 74), (-13, -56), (-14, -55), (6, 14), (7, 13), (5, 15), (-19, -50), (-17, -52), (-18, -51), (-29, -37), (-30, -36), (-24, -44), (-23, -45), (-22, -46), (-20, -48), (-21, -47), (47, 56), (48, 55), (1, 21), (2, 20), (-1, -21), (-2, -20), (29, 37), (30, 36), (1, 40), (2, 39), (3, 38), (-29, -37), (-30, -36), (21, 34), (22, 33), (4, 37), (-21, -34), (-22, -33), (17, 28), (18, 27), (19, 26), (-5, -15), (5, 36), (-19, -26), (-17, -28), (-18, -27), (-6, -14), (-7, -13), (14, 35), (15, 34), (16, 33), (17, 32), (19, 31), (20, 30), (-5, -36), (13, 36), (-4, -37), (12, 37), (-3, -38), (11, 38), (-1, -40), (-2, -39), (8, 40), (9, 39), (-57, -77), (46, 57), (-60, -74), (-58, -76), (-59, -75), (-61, -72), (-62, -71), (-63, -70), (-64, -69), (58, 65), (59, 64), (-46, -57), (6, 70), (7, 69), (5, 71), (4, 72), (1, 75), (2, 74), (-8, -40), (-9, -39), (8, 68), (9, 67), (-47, -56), (-48, -55), (39, 48), (40, 47), (41, 46), (-57, -66), (-58, -65), (-59, -64), (52, 59), (53, 58), (-9, -67), (9, 66), (10, 65), (-41, -46), (-39, -48), (-40, -47), (48, 64), (49, 63), (50, 62), (-9, -66), (-10, -65), (-8, -68), (8, 40), (9, 39), (43, 68), (44, 67), (45, 66), (46, 65), (-50, -62), (-48, -64), (-49, -63), (47, 64), (48, 63)]
    # sequence = 'UAAAAUGAUCACGGUUUCAGCUUUGGACAGGGCGUUCCACUAAACUCCUGGUGACAUAGAUAUAUUGGAUUGCAACUACUCGUCGGUCCGGUUGGCGUUC'
    # s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'
    # s2       = '....(((.(((((((....))).....(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........'


    # # 80 nt, S: -16.60 simple detour, -17.70 better one
    # sequence = "CGGGAGCGAGUAUGGACCGGGUUGAAUUGGUACCUCGGCCGCCCUCUGGGGGUAGCCGACGGGCGGCUUCAGCCGGGCCC"
    # s1 = ".............((.((.(((((((.(.((.(((((((.(((((...))))).))))).)))).).))))))).)).))"
    # s2 = "((....)).....((.((.(((((((...........(((((((((.(.(.....)))).)))))))))))))).))))."

    # # 150 nt, -20.10 kcal/mol direct, indirect -24.50
    # sequence = 'AAAAUAAUGUACCGGACAUUCGCGCACGACCACCAUAUGGCAGAGCAUGUGUCUGUGGACCCACUAUAGCUGGGGCGCUUAACCCCAGAAAAGUAUCUUCGGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAA'
    # s1       = '.............((((((..((...(..(((.....)))..).))..))))))(((((((....(((.((((((.......)))))).....)))....)))))))(((......((.((.........)).))......)))......'
    # s2       = '......(((.((((((...(((....))).(((((((((......))))))...)))((...(((....((((((.......))))))...))).)))))))).)))(((......((.((.........)).))......)))......'

    # sequence = 'AAAAUAAUGUACCGGACAUUCGCGCACGACCACCAUAUGGCAGAGCAUGUGUCUGUGGACCCACUAUAGCUGGGGCGCUUAACCCCAGAAAAGUAUCUUCGGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAA'
    # s1       = '.....(((((.....))))).......(((....(((((......))))))))(((((..(((....((((((((.......)))...............(((...(((.......)))))).........)))))....))).))))).'
    # s2       = '.......((..(((((.....((((.............(((...((.((((...(((((((........((((((.......))))))............)))))))....)))).)).)))........)).))..))).))..))...'

    sequence = 'UAAACUCCUGGUGACAUAGAUAUAUUGGAUUGCAACUACUCGUCGGUCCGGUUGGCGUUCAACCCGCGAUAAGUUGGCGCUUGUCCCACUCCGUAAACCUGUGUCUCUCAGGCGGUUACC'
    s1       = '(((.(((((((.(((((((...(((.(((.((..((....((((((((((((((.....))).))).)))....)))))...))..)).))))))...)))))))..))))).)))))..'
    s2       = '..(((((((((.(((((((...(((.(((((((........).))))))(((((.....)))))...(((((((....))))))).......)))...))))))).).)))).))))...'
    add_moves = [(-23, -95), (-25, -93), (46, 74), (26, 93), (-24, -94)]

    # sequence = 'AGCUUAAUGAUCCACGGGUUCGCUGCACUGCAGCUUUUCACGGACUCUGUGCCGUUUGGAGUGACCAGGUCUCUCCCGCUUCGCGCGAGGUGGGUCCACGCAAGCCGGCACGCGAAAGGG'
    # s1       = '...............(((((((((((...))))).......))))))(((((((((((..((((((...((((...(((...))).))))..))).))).))))).))))))........'
    # s2       = '........((.((...)).))(((((...)))))(((((........(((((((((((..(((....((.(((.(((((.....))).)).)))))))).)))).))))))).)))))..'
    # add_moves = {(22, 41), (68, 97), (12, 119), (69, 96), (-63, -97), (13, 118)}
    # add_moves = {(68, 97), (12, 119), (-61, -99), (69, 96), (-63, -97), (-62, -98), (13, 118)}

    sequence = "CGGGAGCGAGUAUGGACCGGGUUGAAUUGGUACCUCGGCCGCCCUCUGGGGGUAGCCGACGGGCGGCUUCAGCCGGGCCC"
    s1 = ".............((.((.(((((((.(.((.(((((((.(((((...))))).))))).)))).).))))))).)).))"
    s2 = "((....)).....((.((.(((((((...........(((((((((.(.(.....)))).)))))))))))))).))))."
    add_moves = [(57, 62), (56, 63), (55, 64)]

    sequence = "GAAAGUUUUCAAGAAAGUAAUAUAGUGUCCCAGUCAGGGGGUAGUGGCGCGUUCUGGACGCAUAGUAUUC"
    s1 = "(((....))).........(((((((((.((((...(.(.((....)).).).)))))))).)).))).."
    s2 = ".................((.(((.((((.((((..(.(..((....)).).).))))))))))).))..."



    add_moves = []

    # s1, s2 = s2, s1

    section = ()
    # search_width = 64
    search_width = 500
    # search_width = 4200
    Verbose = True
    # Debug = True
    Debug = False

    # add_moves = [(39, 45), (38, 46)]
    # indirect_iterations = 1

    # tabu paper example
    # add_moves = [(1, 43)]
    # add_moves = [(1, 43), (2, 42)]
    # add_moves = [(1, 43), (2, 42), (3, 41)]
    # add_moves = [(1, 43), (2, 42), (3, 41), (4, 40)]
    # indirect_iterations = 1

    # add_moves = []
    indirect_iterations = 1
    indirect_iterations = 2
    # indirect_iterations = 3

    paths = find_path(sequence, s1, s2, indirect_iterations=indirect_iterations, add_moves=add_moves,
                      search_width=search_width, Debug=Debug, Verbose=Verbose)


    search_width = 50
    Verbose = False
    print('orig findpath:')
    pathfinder_result = pathfinder.pathfinder(
        sequence, s1, s2, search_width=search_width, verbose=Verbose)
    print(pathfinder_result.max_en)


    # print (pathfinder.pathfinder(sequence, s2, s1, section=section, search_width=search_width, verbose=Verbose).sE)


def detect_local_minimum(fc, structure):
    # perform gradient walk from sample to determine direct local minimum
    pt = RNA.IntVector(RNA.ptable(structure))
    fc.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
    return RNA.db_from_ptable(list(pt))


