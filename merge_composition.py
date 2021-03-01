#!/usr/bin/env python3
# coding: utf-8

import RNA
import numpy as np

import subprocess
import matplotlib.pyplot as plt
from PIL import Image
from IPython.display import SVG, display
from collections import Counter
from collections import defaultdict

import difflib
import sys
import os
import random
import string
import time
import coloredlogs, logging

import helper
# import pathfinder_python
import pathfinder
# import merge_recursive



def print_d(*args):
    msg = ""
    for m in args:
        msg += str(m)
        msg += " "
    logging.debug(msg)


# coloredlogs.install()

# logging.info("It works!")
# logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')



def merge_check(sequence, s1, s2, Debug=False):

    if Debug:
        coloredlogs.DEFAULT_LOG_FORMAT = '%(levelname)s %(message)s'
        coloredlogs.install(level='DEBUG')

    # logging.info("It works!")
    # logging.debug("This is a log message.")
    # print_d("This is a log message.", "test")
    # logging.error("this is an error message")

    ptables_s1 = RNA.ptable_from_string(s1)
    ptables_s2 = RNA.ptable_from_string(s2)
    ltables_s1 = RNA.loopidx_from_ptable(ptables_s1)
    ltables_s2 = RNA.loopidx_from_ptable(ptables_s2)


    def next_int_loops(min_pos, max_pos):

        '''
        new function to find best interior loops with largest bp_dist
        '''
        print_d("start next int loops", min_pos,max_pos)

        curr_lvl = 0
        diff = 0

        candidates = defaultdict(lambda: [float('inf'),float('inf'),0,0,0])
        c2 = defaultdict(list)
        c3 = []

        for i, (p1, p2) in enumerate(zip(ptables_s1[1:], ptables_s2[1:])):
            # if i==0: continue
            if i < min_pos or i > max_pos: continue

            if p1==0 and p2==0:
                continue
            # if p1==p2 and p1 > i:
            #     curr_lvl += 1
            # elif p1==p2 and p1 < i:
            #     curr_lvl -= 1

            # check which compatible sections have the highest potential for recursion
            if p1==p2 and i<p1:

                j = p1
                # check compatibility
                last_i = i-1
                next_j = j+1

                curr_i = i+1
                next_i = i+2
                last_j = j-1

                # out of bounds
                if next_j>len(s1) or last_i < 0:
                    continue

                # inner/outer section not compatible
                # if ptables_s1[last_i] != ptables_s2[last_i] or ptables_s1[next_j] != ptables_s2[next_j]:        
                #     print_d(last_i,i,j,next_j)
                #     print_d("fail1",ptables_s1[last_i], ptables_s2[last_i] )
                #     print_d("fail2",ptables_s1[next_j], ptables_s2[next_j] )                    
                #     continue

                print_d(curr_i,next_i,last_j,j)
                # print_d("fail1",ptables_s1[curr_i], ptables_s2[curr_i] )
                # print_d("fail2",ptables_s1[i], ptables_s2[i] )  

                if ptables_s1[curr_i] != ptables_s2[curr_i] or ptables_s1[i] != ptables_s2[i]:        
                    print_d(curr_i,next_i,last_j,j)
                    print_d("fail1",ptables_s1[curr_i], ptables_s2[curr_i] )
                    print_d("fail2",ptables_s1[last_j], ptables_s2[last_j] )                    
                    continue
            

                # debatable - extra ( & ) at end/start pos
                # if ptables_s1[last_i] == 0 or ptables_s1[next_j] == 0:  
                #     continue   

                # if start and end are unpaired, both i and j have to be unpaired, 
                # otherwise energies don't add up properly / various errors
            
                last_i = i

                if ptables_s1[last_i] == 0 and ptables_s2[last_i] != 0:
                    continue      
                if ptables_s1[next_j] == 0 and ptables_s2[next_j] != 0: 
                    continue     

                # if ptables_s1[i+1] == 0 or ptables_s1[j] == 0:
                #     continue#
                if ptables_s1[i] == 0 or ptables_s1[j+1] == 0:
                    continue

                # if ptables_s1[i] < ptables_s1[j+1]: # last i,j: )( instead of ()
                if ptables_s1[ptables_s1[i]] != ptables_s1[j+1]: # last i,j: )( instead of ()
                    continue


                outer_s1 = s1[min_pos:i] + "."*(p1-i) + s1[p1:max_pos+1]
                outer_s2 = s2[min_pos:i] + "."*(p2-i) + s2[p2:max_pos+1]



                print_d (outer_s1, s1[i-1], ptables_s1[i], s1[p1], ptables_s1[j+1], "//", ptables_s1[ptables_s1[i]])
                print_d (outer_s2)

                inner_s1 = s1[i:p1]
                inner_s2 = s2[i:p1]

                print_d (inner_s1)
                print_d (inner_s2)

                inner_size = p1-i
                outer_size = max_pos - min_pos - inner_size
                # print (outer_s1)
                # print (outer_s2)
                # print (inner_s1)
                # print (inner_s2)
                bp_dist_inner = RNA.bp_distance(inner_s1, inner_s2)
                bp_dist_outer = RNA.bp_distance(outer_s1, outer_s2)
                # bp_dist = max(bp_dist_outer, bp_dist_inner) - min(bp_dist_outer, bp_dist_inner)

                # optimize = (inner_size/outer_size)*bp_dist_inner
                optimize = (inner_size/outer_size)
                optimize = abs(0.6-(inner_size/outer_size))
                
                # if the step is too small
                # if bp_dist_outer < 3:
                #     continue

                # if min(bp_dist_outer, bp_dist_inner) < 10: # 300_min10
                #     continue
                
                if min(bp_dist_outer, bp_dist_inner) < 3: # 300_r
                    continue

                # this is the standard.. ?
                # if min(bp_dist_outer, bp_dist_inner) < 1:
                #     continue

                # overwrite with better candidate
                print_d ("candidate", i, p1, p2, curr_lvl, diff, "inner size", inner_size, "outer size", outer_size, "opt:", optimize)

                c3.append((i, j, optimize, inner_size, outer_size, bp_dist_inner, bp_dist_outer))

                """
                recursion conditions: inner section > 20 bp
                outer / inner section in the region of 0.25 to 0.75

                maximize inner section size which has at least a 3 bp distance to the outer section

                """

            # print ("---")
        print_d ("found candidates:")

        c3 = sorted(c3, key=lambda item: item[2], reverse=False) # highest opt.
        available = [0 for i in ptables_s1]
        

        for key in c3:
            print_d ("c3", i, j)
            i = key[0]
            j = key[1]

            # all nucleotides between i and j need to be available
            if all(i == 0 for i in available[i+1:j]):
                print_d ('add', i, j)
                available[i:j] = [1] * (j-i)

        mode = 0
        indices = []
        for i in range(len(available)):
            if available[i]==1:
                if mode == 0:
                    start_pos = i
                mode = 1
            else:
                if mode == 1:
                    indices.append([start_pos,i])
                mode = 0
        print_d ("final indices", indices)
        return indices




    def ignore_unpaired_nt(i, j):
        while i<j:
            # if ptables_s1[i] != 0 or ptables_s2[i] != 0:
            if ptables_s1[i] != 0 and ptables_s2[i] != 0:
                # i -= 1
                return i
            i += 1
        # case if no paired nt
        return j


    def ignore_non_aligned(i, j):
        while i<j:
            if ptables_s1[i] == ptables_s2[i] and ptables_s1[i] != 0 and i > ptables_s1[i]:
                # i -= 1
                return i
            i += 1
        # case if no paired nt
        return j
    def bp_dist_section(i,j):
        # bp dist s1 to s2 for a given section
        return RNA.bp_distance(s1[i-1:j], s2[i-1:j])


    def new_exterior_loops(i,j):
        if i==0:
            i = 1
        def exterior_loops_per_s(ptable, ltable, i, j):
            i_s1 = i    
            j_s1 = j
            s1_list = []
            while True:        
                # print ("s", i_s1)
                # ignore unpair nt
                while i_s1<j:
                    if ptable[i_s1] != 0:
                        break
                    i_s1 += 1
                if i_s1 >= j:
                    break
                while j_s1<j:
                    if ptable[j_s1] != 0:
                        break
                    j_s1 -= 1     

                i_s1_start = i_s1
                end_loop = ltable[i_s1]

                i_s1 += 1
                
                while i_s1<j:
                    if ltable[i_s1]==end_loop and ptable[i_s1]!=0:
                        break
                    i_s1 += 1
                # print (i_s1_start,j_s1, i_s1)
                s1_list.append((i_s1_start, i_s1))
                # break
                i_s1 += 1

            return s1_list
        s1_list = exterior_loops_per_s(ptables_s1, ltables_s1, i, j)    
        s2_list = exterior_loops_per_s(ptables_s2, ltables_s2, i, j)
        # print (s1_list)
        # print (s2_list)
        i = 0
        j = 0
        p_min = None
        exterior_loops = []        
        # go over the list of exterior loops of s1 and s2 and see if
        # there are overlapping sections, which have to be merged here

        while (len(s1_list)!=i and len(s2_list)!=j):

            # print ("start", i, j, s1_list[i])

            if p_min == None:
                p_min = min(s1_list[i][0], s2_list[j][0])        

            t1 = s1_list[i]
            t2 = s2_list[j]
            min_t = min(min(t1[0], t2[0]),p_min)
            max_t = max(t1[1], t2[1])
            print_d ('current min/max:', min_t, max_t)

            if i+1 < len(s1_list):
                next_i_min = s1_list[i+1][0]
                print_d ('s1 list:', next_i_min, max_t)
                if next_i_min < max_t:
                    i += 1
                    continue
            if j+1 < len(s2_list):
                next_j_min = s2_list[j+1][0]
                print_d ('s2 list:', next_j_min, max_t)
                if next_j_min < max_t:
                    j += 1
                    continue   
            exterior_loops.append((min_t, max_t))
        
            i += 1
            j += 1
            p_min = None
            # print ("end", i, j)
        # check here for )(  /   ).  /   .(
        #                )(  /   .(  /   ). 
        loop_id = 0
        while len(exterior_loops)>1 and loop_id+1<len(exterior_loops):        
            loop = exterior_loops[loop_id][1]
            next_loop = exterior_loops[loop_id+1][0]
            print_d ("check", loop, next_loop)
            print_d (s1[loop-1], s1[next_loop-1])
            print_d (s2[loop-1], s2[next_loop-1])
            # check for compatibility if end and start for the next loop is next to each other
            # if loop+1 == next_loop and (s1[loop-1] != s2[loop-1] or s1[next_loop-1] != s2[next_loop-1]):
            # if loop+1 == next_loop and (ptables_s1[loop] != ptables_s2[loop] or ptables_s1[next_loop] != ptables_s2[next_loop]):
            if loop+1 <= next_loop and (ptables_s1[loop] != ptables_s2[loop] or ptables_s1[next_loop] != ptables_s2[next_loop]):

                print_d ('merge', loop_id, loop_id+1)
                exterior_loops[loop_id] = (min(exterior_loops[loop_id][0], exterior_loops[loop_id+1][0]),\
                                            max(exterior_loops[loop_id][1], exterior_loops[loop_id+1][1]) )
                exterior_loops.pop(loop_id+1)
                loop_id = 0 
            # overlapping
            if loop>=next_loop:
                print_d ('merge', loop_id, loop_id+1)
                exterior_loops[loop_id] = (min(exterior_loops[loop_id][0], exterior_loops[loop_id+1][0]),\
                                            max(exterior_loops[loop_id][1], exterior_loops[loop_id+1][1]) )
                exterior_loops.pop(loop_id+1)
                loop_id = 0 

            loop_id += 1           
        return exterior_loops




    def recursive_walk(i, j, r_depth=0):
        # collect paths at each recursion level
        sections = []        

        # if r_depth==0:
        #     e = new_exterior_loops(i,j)
        # else:
        #     # e = interior_loops(i,j)
        #     e = next_int_loops(i-1,j-1)
        # print_d ("all ext/int loops:", e)

        e = next_int_loops(i,j)
        for loop in e:        
            ext_i, ext_j = loop
            to_add = [ext_i]

            # if r_depth==0:
            #     return_value = recursive_walk(ext_i, ext_j, r_depth=r_depth+1)
            # else:
            return_value = recursive_walk(ext_i, ext_j-1, r_depth=r_depth+1)
             
            to_add += return_value        
            to_add += [ext_j]
            sections.append(to_add)        

        if r_depth==0:
            sections.insert(0, 0) # start with 0
            sections.append(j) #j-1 # end with total length

        return sections

    return recursive_walk(0,len(s1))






if __name__ == '__main__':


    # Ex 139
    sequence = "AGUUCAGCUACUUACACUGCGGCCUACGGGGGAAUCACUGACCGCGUUAUACCCUAUAUUCGUAUAUGCACGUACAAUACCAACAGUAGGAAUACUUUCUGGGGAGUAGAGUAGAAAUCUUGGCCCCCGACUUGACUGUUAAUCCUCUAC"
    s1 =       ".((...(((...........)))..))(((((..(((......(((.(((((.........))))))))...(((..((((..(((.(((....))).)))..).)))..))).......))))))))(((......))).........."
    s2 =       "...........(((....(((((...((((((..(((......((.((((.(((((.((((.((..((...(((...)))...)).)).))))......))))).)))).))........)))))))))....).)))))))........"
  
    # merge ideas txt
    sequence = "AACGCCGCGCCGGCCAGGCUGCCAGUUGUCUUAUUGAGUUGGCCAUCUACGUAAUACUUUUGUAGCCGACGACUAUAACG"
    s1 =       "(((...(((((.....)).)))..)))...((((...((((((....((((.((....))))))))))))....)))).."
    s2 =       "...((.(((((.....))).))..)).(((.......((((((....((((.........)))))))))))))......."

    #150.90
    #150.4
    sequence = "AGGGGAAAGUUAUCGAACAGACCUUGUUUGGCUUAUGGCACUUGGGUAGAUUUUCUCUACCGUUGUCAUGAAAACGCACCGCAAGAGAAGUCACGAACCUGCGCCAGCGCACUACAGCUUGACUUUGGGAACUGAGUGGCGUAGUGACGG"
    s1 =       "...((...(((.(((((((.....))))))).((((((((....((((((.....))))))..)))))))).)))...)).........(((((.....(((((((((........)))..((((.(....).))))))))))))))).."
    s2 =       "...((...((..(((.......(((((...((((((((((..((.(((((.....))))))).))))))))....))...)))))........)))....)).))..((((((((.(((..((((........))))))))))))).))."
    # sections: [0, [39, 64], [124, 134], 149]

    # 60.7
    sequence = "UUGACCAUUGACUACGGGAAAUUUUCUGUUUCGUUCUCCAGGGAGGCCAAGGAGCUGUCC"
    s1 =       "((((...))))..((((((....))))))..((..((((..(.....)..)))).))..."
    s2 =       "...........(....)(((((.....)))))...((((.((....))..))))......"

    # # 60.29
    # sequence = "GAGGUGUAGUCCAGCCCUUAUUUGUUACCCGGCCGAUCCGAGAUAUGCACGUACGUCCAA"
    # s1       = "((.(((((....(((........)))..(((.......)).)...))))).).)......"
    # s2       = "((((.((......))))))..........(((.....)))...................."


    # # 60.49
    # sequence = "GAAGGCGGUUGGGGUGCAUACCUGCCAUACCAAUUCGUAUAUUUCGUCGUCUCCUCCUGC"
    # s1       = "((((((((((((..((((....))).)..))))).)))...))))..............."
    # s2       = "(.((((((((((...((......))....))))).((.......)).)))))...)...."


    sequence = "GAAGUGUGGCUGGCAAGGAGAAUUAUGUGUGAAAAUUUGUCGGUAGAUAGGCAGUGGUGGCGAAGGGAGGGGGAAACGAUUUUGCCUCCGACGUCCAUCAUCGCCAGACAGGACGGUCUCCCUUCCUACAGGUCUCUGGCACAUAUCCUC"
    s1       = "....(((.....))).(((....(((((((.(.....((((....))))(((.((((((((....(((((.(((......))).)))))...).))))))).)))((((((((.((....)).))))....)))).).)))))))))).."
    s2       = "((.(((((((.((.(((((((....(.(((.......((((....))))(((.(((((((((...(((((.(((......))).)))))..)).))))))).)))..))).)....)))).))))).....)))....))))...))..."

    sequence = "GGAUAUUUCUUGUUGGCGCUCGGGCCGUCACUCUCCUCCCAACGAAACCCCAGGAGAGACAUCACAUAAGCAAACCUUUUGAUUUGAUGUAACCGUGGAGAAAACAAGUUCCUGUUACUUGGACACGUCUUUAGAAAAAACAGGAACGGU"
    s1       = "......((((.((((((......)))....(((((((..............)))))))((((((((.............))...)))))))))...))))...((..(((((((((....(((....)))........))))))))).))"
    s2       = "........(((((((((......)))....(((((((..............)))))))((((((.((((.........)).)).)))))).....((((((...(((((.......)))))......)))))).....))))))......"

    sequence = "AACGGGUGGGUACUCCCUGGUAAAGCCCGAGUCGAGACAUUGUCAUAUGUAUGAGAUUCCUUUGUUGUUGGUCGGCUGGG";
    s1       = "..((((((((....))).......)))))((((((.....(.((((....)))).)...((........))))))))...";
    s2       = "...((((...((((....))))..)))).(((((((((....((..........)).......))).....))))))...";


    Debug = True
    # Debug = False

    sections = merge_check(sequence, s1, s2, Debug=Debug)

    print ("input")
    print (sequence)
    print (s1)
    print (s2)

    fc = RNA.fold_compound(sequence)
    pt1 = RNA.ptable_from_string(s1)
    pt2 = RNA.ptable_from_string(s2)

    helper.print_tables(s1,s2,pt1,pt2)
    # res = merge_recursive.recursive_merge(sequence, s1, s2, sections=sections, search_width=500, Debug=False, Verbose=True, new=True, plot_graph=False)


    print (sections)