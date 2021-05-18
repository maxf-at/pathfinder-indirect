#!/usr/bin/env python3
# coding: utf-8

# generic Python libraries
import numpy as np


from collections import defaultdict
from dataclasses import dataclass

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
    moves_indirect: list

    # energies:     list
    s_hash = 0
    opt = 0
    add_moves = []


class Fp_class():
    def __init__(self, sequence, s1, s2, search_width, add_moves=[]):

        self.fc = RNA.fold_compound(sequence)
        self.sequence = sequence
        self.s1 = s1
        self.s2 = s2
        self.search_width = search_width

        self.moves = []

        self.h_list = []
        self.bp_dist = 0
        self.add_moves = add_moves

        self.init()

    def int_hash_64(self, x) -> np.uint64:

        z = np.uint64(x) + 0x9e3779b97f4a7c15
        # z = (x += 0x9e3779b97f4a7c15)
        z = (z ^ (z >> np.uint64(30))) * 0xbf58476d1ce4e5b9
        z = (z ^ (z >> np.uint64(27))) * 0x94d049bb133111eb
        return z ^ (z >> np.uint64(31))

    def init(self):
        np.set_printoptions(suppress=True)  # prevent numpy exponential

        max_energy = float("inf")

        current_width = self.search_width  # for the time being, only 1 pass

        s1 = self.s1
        s2 = self.s2

        print(s1)
        print(s2)

        pt1 = list(RNA.ptable_from_string(s1))
        pt2 = list(RNA.ptable_from_string(s2))
        self.bp_dist = RNA.bp_distance(s1, s2)

        move_list = []
        h_list = []

        end_hash = np.uint64(0)  # allow controlled overflow
        bp_dist = 0

        for i in range(1, pt1[0]):
            if (pt1[i] != pt2[i]):
                if (i < pt1[i]):
                    move_list.append((-i, -pt1[i]))

                    move_hash = self.int_hash_64(-i*pt1[i])
                    h_list.append(move_hash)
                    # h_list = np.append(h_list, np.uint64(move_hash))

                    end_hash += np.uint64(move_hash)

                    bp_dist += 1

                if (i < pt2[i]):
                    move_list.append((i, pt2[i]))

                    move_hash = self.int_hash_64(i*pt2[i])
                    h_list.append(move_hash)

                    end_hash += np.uint64(move_hash)

                    bp_dist += 1

        # print (move_list)
        # print (np.array(h_list).astype(np.uint64))
        # print (h_list)

        # print (end_hash)
        # print (np.sum(h_list))

        self.moves = move_list
        self.h_list = h_list

        path_generator = self.find_path_once(
            s1, s2, max_energy, current_width, mode=True, sort_min=False, Debug=Debug, Verbose=Verbose)

        for path in path_generator:

            # current_moves =

            current_moves = [0] * bp_dist

            for a, b in enumerate(path.moves):
                current_moves[b] = self.moves[a]

            print(path)
            print(current_moves)
            print(self.print_moves(sequence, s1, s2, current_moves, path.moves_indirect))
            break

        # # iterate over all currently found paths. somepaths might not be optimal
        # for path in path_generator:

    def find_moves(self, temp_pt, current_moves, temp_loop):
        """
        generator function, yields possible structures 1 move away
        from the original structure by finding fitting i and j with
        RNA pair and loop tables

        """
        # loop table
        # temp_loop = RNA.loopidx_from_ptable(temp_pt)
        c = False

        for a in range(len(current_moves)):

            if current_moves[a] > 0:
                continue

            i = self.moves[a][0]
            j = self.moves[a][1]

            if (j < 0 and temp_pt[-i] == -j):

                # c = alt_pt[-i] == -j # compatibility with alternative ptable
                yield a, i, j

            else:
                if (j > 0 and (temp_loop[i] == temp_loop[j]) and
                        (temp_pt[i] == 0) and (temp_pt[j] == 0)):

                    # c = (alt_pt[i] == 0) and (alt_pt[j] == 0)
                    yield a, i, j

    def find_path_once(self, s1, s2, max_energy, width, mode=True, sort_min=False, Debug=False, Verbose=False):
        """
        main findpath algorithm (bounded BFS)
        """

        #
        add_moves = self.add_moves

        if Debug:
            coloredlogs.DEFAULT_LOG_FORMAT = '%(levelname)s %(message)s'
            coloredlogs.install(level='DEBUG')

        bp_dist = self.bp_dist

        pt1 = list(RNA.ptable_from_string(s1))
        pt2 = list(RNA.ptable_from_string(s2))
        fc = RNA.fold_compound(self.sequence)

        e1 = fc.eval_structure_pt(pt1)

        current_bp = 0
        current_bp_end = bp_dist

        empty_moves = [-1] * bp_dist
        empty_moves_i = [set()] * bp_dist

        init_intermediate = Intermediate(p_table=pt1.copy(), mode=mode, saddle_e=e1, current_e=e1, moves=empty_moves.copy(), moves_indirect=empty_moves_i.copy())

        # paths start with 1 initial intermediate
        paths = [init_intermediate]

        # dont stop at current_bp_end, consider potential indirect moves
        while (current_bp != current_bp_end):

            # collect all new paths here (next iteration)
            collect_paths = []

            # if current_bp == 2:
            #     break

            for current_path in paths:

                current_p_table = current_path.p_table
                current_s = RNA.db_from_ptable(current_p_table)

                current_e = current_path.current_e
                current_saddle_e = current_path.saddle_e
                current_s_hash = current_path.s_hash
                current_moves = current_path.moves
                current_moves_indirect = current_path.moves_indirect

                current_indirect_state = set()
                for i, j in add_moves:
                    if i < 0:
                        continue
                    if current_p_table[i] == j:
                        current_indirect_state.add((i, j))

                current_moves_indirect[current_bp] = current_indirect_state

                if current_bp < 5:
                    print("d:", current_bp, current_s, current_saddle_e, current_indirect_state)

                temp_loop = RNA.loopidx_from_ptable(current_p_table)

                # "try_moves"
                for a, i, j, in self.find_moves(current_p_table, current_moves, temp_loop):
                    # print ("t", a, i, j)
                    # next energy calculations
                    next_e = fc.eval_move_pt(current_p_table, i, j) + current_e

                    # print ("en done")
                    next_p_table = current_p_table.copy()
                    if i < 0:
                        next_p_table[-i] = 0
                        next_p_table[-j] = 0
                    else:
                        next_p_table[i] = j
                        next_p_table[j] = i
                    # print ("d:", current_bp, RNA.db_from_ptable(next_p_table), next_e)
                    # next_e = min(next_e, next_e_alt)
                    # next saddle energy
                    next_s = max(current_saddle_e, next_e)
                    # if this move is valid... append to list
                    if next_s <= max_energy:
                        # print ("en done")
                        next_moves = current_moves.copy()
                        next_moves[a] = current_bp
                        new_intermediate = Intermediate(p_table=next_p_table, mode=mode, saddle_e=next_s, current_e=next_e,
                                                        moves=next_moves, moves_indirect=current_moves_indirect.copy())
                        collect_paths.append(new_intermediate)

                # second try moves iteration, alternatives

                # try to adjust the current structure

                # check if we should remove or add stabilizers
                pos = 0
                done = 0
                alt = False

                candidate_extra_move = (0, 0)
                candidate_extra_move_en = 9999

                for i, j in add_moves:
                    if i < 0 and current_p_table[-i] == -j and current_bp > 23:
                        alt = True

                        test_en = fc.eval_move_pt(current_p_table, i, j) + current_e

                        print("remove", i, j, test_en)

                        if test_en < candidate_extra_move_en:
                            candidate_extra_move_en = test_en
                            candidate_extra_move = (i, j, test_en)

                        # current_p_table[-i] = 0
                        # current_p_table[-j] = 0
                        # current_indirect_state.remove((-i, -j))

                    # only add if we're on the same loop
                    if i > 0 and current_p_table[i] == 0 and current_p_table[j] == 0\
                         and temp_loop[i] == temp_loop[j] and current_bp < 25:

                        alt = True

                        test_en = fc.eval_move_pt(current_p_table, i, j) + current_e
                        if test_en < candidate_extra_move_en:
                            candidate_extra_move_en = test_en
                            candidate_extra_move = (i, j)

                        # current_p_table[i] = j
                        # current_p_table[j] = i
                        # current_indirect_state.add((i, j))

                        # print ("add",i,j, RNA.db_from_ptable(current_p_table))

                current_moves_indirect[current_bp] = current_indirect_state

                if alt == False:
                    continue

                # adjust ptable with indirect move

                i = candidate_extra_move[0]
                j = candidate_extra_move[1]

                if i > 0:
                    current_p_table[i] = j
                    current_p_table[j] = i
                    current_indirect_state.add((i, j))
                else:
                    current_p_table[-i] = 0
                    current_p_table[-j] = 0
                    current_indirect_state.remove((-i, -j))

                temp_loop = RNA.loopidx_from_ptable(current_p_table)
                current_s = RNA.db_from_ptable(current_p_table)

                # current_e = fc.eval_structure_pt(current_p_table)
                current_e = candidate_extra_move_en
                current_saddle_e = max(current_e, current_saddle_e)

                # if current_bp <5:
                print("best c:", candidate_extra_move, candidate_extra_move_en)
                print("a:", current_bp, current_s, current_e, current_indirect_state)

                # "try_moves"
                for a, i, j, in self.find_moves(current_p_table, current_moves, temp_loop):
                    # print ("t", a, i, j)
                    # next energy calculations
                    next_e = fc.eval_move_pt(current_p_table, i, j) + current_e

                    # print ("en done")
                    next_p_table = current_p_table.copy()
                    if i < 0:
                        next_p_table[-i] = 0
                        next_p_table[-j] = 0
                    else:
                        next_p_table[i] = j
                        next_p_table[j] = i
                    # print ("d:", current_bp, RNA.db_from_ptable(next_p_table), next_e)
                    # next_e = min(next_e, next_e_alt)
                    # next saddle energy
                    next_s = max(current_saddle_e, next_e)
                    # if this move is valid... append to list
                    if next_s <= max_energy:
                        # print ("en done")
                        next_moves = current_moves.copy()
                        next_moves[a] = current_bp
                        new_intermediate = Intermediate(p_table=next_p_table, mode=mode, saddle_e=next_s, current_e=next_e,
                                                        moves=next_moves, moves_indirect=current_moves_indirect.copy())
                        collect_paths.append(new_intermediate)

            # print (current_bp, "sort")

            # first sorting step
            collect_paths.sort(key=lambda x: (x.p_table, x.saddle_e))

            last_ptable = []
            # last_ptable = collect_paths[-1][0].p_table
            # last_ptable = collect_paths[0][-1].p_table
            print_d("sort done", last_ptable, init_intermediate.p_table)

            # remove duplicates ptables
            if current_bp+1 != current_bp_end:
                for i in range(len(collect_paths)):
                    c_ptable = collect_paths[i].p_table
                    if c_ptable == last_ptable:
                        # set saddle energy high
                        collect_paths[i].saddle_e = 999
                    else:
                        last_ptable = c_ptable
                i = 0
                while True:
                    if i > len(collect_paths)-1:
                        break
                    if collect_paths[i].saddle_e == 999:
                        collect_paths.pop(i)
                        continue
                    i += 1

            # second sorting step
            collect_paths.sort(key=lambda x: (x.saddle_e, x.current_e))

            # discard paths
            collect_paths = collect_paths[:width]

            # return valid paths if we're past the minimum bp_dist
            # if current_bp >= current_bp_end-1:
            #     for i in range(len(collect_paths)):
            #         if collect_paths[i][-1].p_table == list(end_p_table):
            #             yield collect_paths[i][-1]

            # next iteration
            paths = collect_paths
            current_bp += 1

        # return remaining paths
        if paths:
            for path in paths:
                if path:
                    yield path

    def print_moves(self, sequence, s1, s2, moves, indirect_moves, move_color='\033[93m', Verbose=True, ignore_moves=set()):
        """
        pretty print a folding path with colour coding

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

        if Verbose:
            print(sequence)
        s = s1
        pt = list(RNA.ptable_from_string(s))

        move_color = '\033[92m'

        fc = RNA.fold_compound(sequence)
        e1 = en = round(fc.eval_structure(s), 2)
        max_en = float("-inf")

        output_rows = []
        moves_i_j = [(x[0], x[1]) for x in moves]

        last_indirect_state = set()

        # preprocessing - generate strings & energies if required
        for a, move in enumerate(moves):
            en = False


            # if indirect_moves[a] != set():

            new_move = indirect_moves[a] - last_indirect_state
            new_move2 = list(last_indirect_state - indirect_moves[a])
            if new_move != set():
                new_move = list(new_move)[0]
            if new_move2 != []:
                # undo indirect moves
                new_move = (-new_move2[0][0], -new_move2[0][1])
            last_indirect_state = indirect_moves[a]





            en2 = -9999
            if new_move != set():
                i = new_move[0]
                j = new_move[1]
                if i > 0:  # add bp
                    pt[i] = j
                    pt[j] = i
                if j < 0:  # del bp
                    pt[-i] = 0
                    pt[-j] = 0
                s = RNA.db_from_ptable(pt)
                en2 = round(fc.eval_structure(s), 2)


            if len(move) == 2:
                i, j = move
            if len(move) == 3:
                i, j, en = move

            if (i, j) in ignore_moves:
                continue

            lp = RNA.loopidx_from_ptable(pt)

            # print (i,j, lp[i], lp[j])

            if i > 0:  # add bp
                pt[i] = j
                pt[j] = i

            if j < 0:  # del bp
                pt[-i] = 0
                pt[-j] = 0



            s = RNA.db_from_ptable(pt)
            print (s)
            if not en:
                en = round(fc.eval_structure(s), 2)

            en = max(en, en2)

            e2 = en
            if en > max_en:
                max_en = en
            output_rows.append((s, i, j, en))

        for s, i, j, en in output_rows:

            # print initial row with move (0,0)
            if i == 0:
                info = f'{move_color}[{i:4}, {j:4} ]{c.ENDC} {en:6.2f}'
                if Verbose:
                    print(f"{s} {info}")
                continue

            pos_i = abs(i)-1
            pos_j = abs(j)-1

            # if a move has an identical inverse copy (1,2) <-> (-1,-2)
            # it is automatically an indirect move - these are colored in red
            if (-i, -j) in moves_i_j:  # indirect move
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

            if Verbose:
                print(f"{colored_s} {info}")

        barrier = max_en - e1
        if Verbose:
            print(
                f"S: {max_en:6.2f} kcal/mol | B: {barrier:6.2f} kcal/mol | E[start]:{e1:6.2f} E[end]:{e2:6.2f}")

        if ignore_moves != set():
            return max_en, s  # last structure

        return max_en


if __name__ == '__main__':

    # various random sequences

    # # 60.1 inner
    # sequence = "GCCAACAAACCGUGAUGGGCUAUGUUUGAUUCAUCCUAUUUAUGUUUUUCGAGAUGCGCG"
    # s1       = "...........(((..(((..((.....))...)))....)))................."
    # s2       = "...........((((((((..(((.......)))))))).)))................."

    # ~~~

    # Tabu paper example
    # sequence = "CGCGACGGCUACGCGACGGCAAUGCCGUUGCGAAGCCGUCGCGAUC"
    # s1 = "(((((((((..............))))))))).............."
    # s2 = "...........(((((((((..............)))))))))..."

    # rna2dfold example
    sequence = "GGGCGCGGUUCGCCCUCCGCUAAAUGCGGAAGAUAAAUUGUGUCU"
    s1 = "(((((.....)))))(((((.....)))))(((((.....)))))"
    s2 = "((((((((((.....(((((.....))))).....))))))))))"
    add_moves = set([(6, 12), (-6, -12), (7, 11), (-7, -11), (5, 13), (-5, -13)])

    sequence = "CGGGAGCGAGUAUGGACCGGGUUGAAUUGGUACCUCGGCCGCCCUCUGGGGGUAGCCGACGGGCGGCUUCAGCCGGGCCC"
    s1 = ".............((.((.(((((((.(.((.(((((((.(((((...))))).))))).)))).).))))))).)).))"
    s2 = "((....)).....((.((.(((((((...........(((((((((.(.(.....)))).)))))))))))))).))))."
    add_moves = [(57, 62), (-57, -62), (56, 63), (-56, -63), (55, 64), (-55, -64)]

    # s1, s2 = s2, s1

    section = ()
    search_width = 5000
    Verbose = True
    # Debug = True
    Debug = False

    # add_moves = [(39, 45), (38, 46)]
    # indirect_iterations = 1

    # add_moves = []
    # indirect_iterations = 2
    # indirect_iterations = 2

    fp = Fp_class(sequence, s1, s2, search_width, add_moves)

    # print (se)
    print('orig findpath:')

    pathfinder_result = pathfinder.pathfinder(
        sequence, s1, s2, search_width=search_width, verbose=Verbose)

    print(pathfinder_result.max_en)
    # print (pathfinder.pathfinder(sequence, s2, s1, section=section, search_width=search_width, verbose=Verbose).sE)
