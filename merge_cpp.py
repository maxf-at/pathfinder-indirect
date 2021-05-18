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

from helper import path_class

import pathfinder_cpp
import pathfinder
import merge_composition

import networkx as nx

def print_d(*args):
    msg = ""
    for m in args:
        msg += str(m)
        msg += " "
    # logging.debug(msg)

def p_to_s(tb):
    """ptable --> string"""
    s = ""
    for i, p in enumerate(tb):
        if i==0: continue
        if p == 0:
            s += "."
        elif p>i:
            s += "("
        else:
            s += ")"
    return s


process_label = lambda x: float(x.split("/")[1])

def permutations(G_i, G_j, current_node_i, current_node_j, section_i, section_j, mode):
    """
    returns possible future path permutations from 2 input digraphs
    as list of tuples, e.g.
    [(1, 0), (2, 0), (3, 0), (3, 1), (3, 2), (3, 3)]
    """

    if mode: # fwd pass
        e_2 = process_label(G_j.nodes[current_node_j]['label'])
        # for edge_1 in G_i[current_node_i]:
        for _, edge_1 in G_i.out_edges(current_node_i):
            e_1 = process_label(G_i.nodes[edge_1]['label'])
            move = G_i[current_node_i][edge_1]['label']
            yield (edge_1, current_node_j, e_1+e_2, move)

        e_1 = process_label(G_i.nodes[current_node_i]['label'])
        # for edge_2 in G_j[current_node_j]:
        for _, edge_2 in G_j.out_edges(current_node_j):
            e_2 = process_label(G_j.nodes[edge_2]['label'])
            move = G_j[current_node_j][edge_2]['label']
            yield (current_node_i, edge_2, e_1+e_2, move)

    else: # bwd pass
        # print (current_node_i, current_node_j)
        e_2 = process_label(G_j.nodes[current_node_j]['label'])
        # for edge_1 in G_i[current_node_i]:
        for edge_1, _ in G_i.in_edges(current_node_i):
            # print (current_node_i, edge_1)
            e_1 = process_label(G_i.nodes[edge_1]['label'])
            move = G_i[edge_1][current_node_i]['label']
            yield (edge_1, current_node_j, e_1+e_2, move)

        e_1 = process_label(G_i.nodes[current_node_i]['label'])
        # for edge_2 in G_j[current_node_j]:
        for edge_2, _ in G_j.in_edges(current_node_j):
            e_2 = process_label(G_j.nodes[edge_2]['label'])
            move = G_j[edge_2][current_node_j]['label']
            yield (current_node_i, edge_2, e_1+e_2, move)



def merge_findpath(sequence, s1, s2, paths_i, paths_j, section_i, section_j, merge_search_width=1, mid=False, Debug=False, Verbose=False):

    """
    Merge 2 sections with findpath lookahead
    either inner/outer or adjacent sections
    """

    if Debug:
        coloredlogs.DEFAULT_LOG_FORMAT = '%(levelname)s %(message)s'
        coloredlogs.install(level='DEBUG')

    start = time.time()

    # paths = []
    # for section in sections:
    #     path, max_energy = pathfinder_python.find_path(sequence, s1, s2, section, Verbose=True)
    #     paths.append(path)
    #     print (path)
    print_d(sequence)
    print_d(s1)
    print_d(s2)
    print_d (section_i, section_j)

    G_i = paths_i.G
    G_j = paths_j.G


    # print ("debug:")
    # print (G_i.number_of_nodes())
    # print (nx.adjacency_matrix(G_i).todense())
    # for i_n in G_i.nodes():
    #     for j_n in G_j.nodes():
    #         print (i_n[:-1])
    #         break



    mode = True # fwd pass, false is bwd pass

    # plot_graph(G_i, "i_graph.png")
    # plot_graph(G_j, "j_graph.png")

    if Verbose:
        print ("Start merging, number of nodes:", G_i.number_of_nodes(), G_j.number_of_nodes(), paths_i.dist)




    # iterations = [1,2,4,8,16,32,64,128,256] # maxL
    # iterations = [100, 100, 100]
    iterations = [2, 100, 100]
    # iterations = [2, 50, 50]
    # iterations = [2, 20,20]
    # merge_search_width = None

    max_en = float('inf')

    return_paths = path_class()
    return_paths.s1 = s1
    return_paths.s2 = s2

    results_max_en = []




    for size in iterations:

        max_en = float('inf')
    #   

        print_d ("mode:", mode)

            # current indices (digraph?)
        current_i = 0
        current_j = 0
        i_end = paths_i.dist
        j_end = paths_j.dist

        fc = RNA.fold_compound(sequence)
        
        current_bp_dist = 0

        if mode:
            current_ptable = list(RNA.ptable_from_string(s1))
            current_node_i = list(G_i.nodes)[0]
            current_node_j = list(G_j.nodes)[0]
 
        else:
            current_ptable = list(RNA.ptable_from_string(s2))
            current_node_i = paths_i.s2+str(paths_i.dist)
            current_node_j = paths_j.s2+str(paths_j.dist)

        current_e = fc.eval_structure_pt(current_ptable)/100
        bp_dist_end = i_end+j_end

        if size>10 and merge_search_width:
            size = int(bp_dist_end*merge_search_width)
            if size<11:
                size = 15
            # print ("setting merge search width to", size)

        # print ("current merge search width:", size)

        paths = [[[current_ptable,current_e,current_e,0,0,current_node_i, current_node_j]]]

        all_paths = [] # at end for all fwd and bwd paths

        # print ('start', mode, bp_dist_end, 'current e', current_e, current_node_i[-1], current_node_j[-1])

        while (current_bp_dist != bp_dist_end):

            # if Verbose:
            #     print ("iteration", size, current_bp_dist)


            # print_d (current_bp_dist, mode, bp_dist_end)
            collect_paths = []

            for current_path in paths:
                
                current_ptable,last_s,last_e,_,_,current_node_i, current_node_j = current_path[-1]
                print_d("start perm", current_node_i, current_node_j)


                for (next_node_i, next_node_j, en, move) in permutations(G_i, G_j, current_node_i, current_node_j, section_i, section_j, mode):
                    print_d ("d =", current_bp_dist, next_node_i[-20:], next_node_j[-20:], en, move)
                    

                    new_ptable = current_ptable.copy()
                    
                    l_table = RNA.loopidx_from_ptable(new_ptable)

                    # print_d(RNA.db_from_ptable(new_ptable))

                    # if True:
                    i = move[0]
                    j = move[1]
                    if i==0:
                        continue
                    if mode: # fwd pass: delete bp if neg.
                        if i<0:
                            if (new_ptable[-i] != -j):
                                continue
                            print_d("a")
                            new_ptable[-i] = 0
                            new_ptable[-j] = 0
                        else:
                            if (new_ptable[i] != 0 or new_ptable[j] != 0):
                                continue
                            if (l_table[i] != l_table[j]):
                                continue
                            print_d("b")
                            new_ptable[i] = j
                            new_ptable[j] = i  
                    if not mode: # bwd pass: add bp if neg.
                        if i<0:
                            if (new_ptable[-i] != 0 or new_ptable[-j] != 0):
                                continue
                            if (l_table[-i] != l_table[-j]):
                                continue
                            print_d("c", new_ptable[-i], new_ptable[-j], i, j)
                            new_ptable[-i] = -j
                            new_ptable[-j] = -i
                        else:
                            if (new_ptable[i] != j):
                                continue
                            print_d("d")
                            new_ptable[i] = 0
                            new_ptable[j] = 0

                    # print_d(RNA.db_from_ptable(new_ptable))

                    current_e = fc.eval_structure_pt(new_ptable)/100
                    current_s = max(last_s, current_e)

                    if current_e <= max_en:
                        new_path = current_path.copy() + [[new_ptable,current_s, current_e,move[0],move[1],next_node_i, next_node_j]]
                        print_d("append")
                        collect_paths.append(new_path)


            collect_paths.sort(key=lambda x: (x[-1][0],x[-1][1])) # sort last elem which was just added

            last_ptable = collect_paths[-1][0] # should be [0][-1]??
            # remove duplicates ptables
            for i in range(len(collect_paths)):
                c_ptable = collect_paths[i][-1][0]
                if c_ptable == last_ptable:
                    collect_paths[i][-1][1] = 999 # set saddle energy high
                else:
                    last_ptable = c_ptable

            collect_paths.sort(key=lambda x: (x[-1][1],x[-1][2])) # sort by saddle en, current en

            # for path in collect_paths[:5]:
            #     print (p_to_s(path[-1][0]), path[-1][1:3]) # last elem without ptable

            collect_paths = collect_paths[:size]

            paths = collect_paths

            current_bp_dist += 1

        if size>9 and mode: # last round
            for current_ptable,current_s,current_e,i,j,current_node_i, current_node_j in paths[0]:
                if Verbose:
                    print (f'{p_to_s(current_ptable)} {current_e:5} {current_s:5} ({i:3}/{j:3})')
            return_paths.max_en = current_s            
            return_paths.dist = paths_i.dist + paths_j.dist
            
            all_paths += paths

        # if size == iterations[-1] and not mode: # last round bwd
        if size >9 and not mode: # last round bwd
            last_i, last_j = 0, 0
            for current_ptable,current_s,current_e,i,j,current_node_i, current_node_j in reversed(paths[0]):
                # if Verbose:
                    # print (f'{p_to_s(current_ptable)} {current_e:5} {current_s:5} ({last_i:3}/{last_j:3})')
                last_i = i
                last_j = j
            return_paths.max_en = current_s            
            return_paths.dist = paths_i.dist + paths_j.dist

            bwd_paths = []
            for path in paths:
                current_path = []
                last_i, last_j = 0, 0
                for current_ptable,current_s,current_e,i,j,current_node_i, current_node_j in reversed(path):
                    current_path.append([current_ptable,current_s,current_e,last_i,last_j,current_node_i, current_node_j])
                    last_i, last_j = i, j
                bwd_paths.append(current_path)

        max_en = paths[0][-1][1]
        results_max_en.append(max_en)
        # if Verbose:
        #     print ("Merge Step Result, S:", max_en, "merge runtime:", round(time.time()-start,2), "s")
        
        # switch from fwd to bwd and other way round
        mode = not mode

    if Verbose:
        print ("Merge Step Result, S:", max_en, "merge runtime:", round(time.time()-start,2), "s")


    return_paths.runtime_merge = time.time()-start



    # extract (i,j) moves for display
    


    all_paths += bwd_paths

    all_paths.sort(key=lambda x: (x[-1][1],x[-1][2])) # sort by saddle en, current en

    return_paths.paths = [(i[3],i[4]) for i in all_paths[0]]

    # print ("fwd:")
    # for move in paths[0]:
    #     print (move)
    # print ("bwd:")
    # for move in bwd_paths[0]:
    #     print (move)


    return_paths.max_en = min(results_max_en)

    G = nx.DiGraph()
    # G = nx.Graph()

    current_ptable = list(RNA.ptable_from_string(s1))
    e = fc.eval_structure_pt(current_ptable)/100
    last_node = p_to_s(current_ptable)+str(0)
    first_node = last_node
    G.add_node(last_node)
    G.nodes[last_node]['color'] = "g"
    G.nodes[last_node]['label'] = str(0)+" / "+str(e)

    # process best paths into nodes and edges...
    for path in all_paths:

        last_node = first_node

        for l, move in enumerate(path[1:]):
            current_ptable = move[0]
            e = move[1]
            s = move[2]
            i = move[3]
            j = move[4]
            # print (i,j)
            current_node = p_to_s(current_ptable)+str(l+1)
            e = round(fc.eval_structure_pt(current_ptable)/100,2)
            G.add_node(current_node)
            if l==return_paths.dist-1:
                G.nodes[current_node]['color'] = "g" # end node
            else:
                G.nodes[current_node]['color'] = "r"
            G.nodes[current_node]['label'] = str(l+1)+" / "+str(e)
            G.add_edge(last_node, current_node)
            G[last_node][current_node]['label'] = (i,j)
            G[last_node][current_node]['weight'] = 1

            last_node = current_node

    return_paths.G = G

    return_paths.bp_dist = [paths_i.bp_dist, paths_j.bp_dist]

    return_paths.runtime_ext = time.time()-start


    if not Debug:
        return return_paths


def plot_graph(G, filename, x=20, y=12):

    G = G.to_undirected()
    
    x = 12
    y = 8

    # visualization below
    # plt.figure(figsize=(12, 8), dpi=250, facecolor='w', edgecolor='k')
    plt.figure(figsize=(x, y), dpi=250, facecolor='w', edgecolor='k')
    pos = nx.spring_layout(G, dim=2, iterations=50)
    edge_labels = nx.get_edge_attributes(G, "label")
    node_labels = nx.get_node_attributes(G, "label")

    # label offset
    pos_attrs = {}
    for node, coords in pos.items():
        pos_attrs[node] = (coords[0], coords[1] - 0.07)

    formatted_edge_labels = {(elem[0],elem[1]):edge_labels[elem] for elem in edge_labels} # use this to modify the tuple keyed dict if it has > 2 elements, else ignore

    colors = [node[1]['color'] for node in G.nodes(data=True)]
    nx.draw(G, pos, node_color=colors)
    nx.draw_networkx_labels(G,pos_attrs,node_labels,font_size=15,font_color='b')
    nx.draw_networkx_edge_labels(G,pos,edge_labels=formatted_edge_labels,font_size=11, font_color='red')

    # plt.show()
    plt.savefig(filename)



if __name__ == '__main__':


    # really good example, where old findpath fails...
    # sequence = 'AUGUUCUCAGUUGUUGUUGGAGCAGUCCCGUAGGAGGUCCCCAUUGGGCCCUUUUCUAGGGCUGCGAAGUCGGGCUCCAG'
    # s1 =       '.((....)).......((((((((((((....((((((((.....))))).)))....))))))((....))..))))))'
    # s2 =       '...........((..(((.((((((.(((.(((((((.(((....))).)))..)))))))))))....)).)))..)).'
    
    # ex 50
    # sequence = 'UCAUCCAAAGAACUCGUCAUGGGCGUCAACUUCAGUCGGUAAAGCUCACAGGUGACUUUGAUUUAAGCAAUGAGAUCUUGGACCUCUGCUAGAAGCUAUGACCUUGUCUGCCAAGCUUUUACUCAAAAAGACUCGUAUUCUCAAUCGUGA'
    # s1 =       '...(((((.((.(((((.......(((((..(((.(..((.......)).).)))..))))).......))))).))))))).....(((...)))(((((.......(((..(((((((.....))))).)).))).......))))).'
    # s2 =       '...(((.((((.(((((.......(((.((((..((.((.....)).))))))))).............))))).))))))).......((((((((..(((...))).....)))))))).......((.........)).........'
    # section_1 = [3, 5, 80, 150]
    # section_2 = [5, 80]


    # Ex 75 pos
    # section_1 = [10, 29, 131, 149]
    # section_2 = [29, 131]
    # sequence = "CGCCACGAACACUAUGCGCAACCUACAAGGUCUUCCCCCUCCGGCUGGUAGCAUACACAUGUACACAUCUUCGUUCAAAAUGGCACGGAACGGCCAAGAGAAGCACCGCUACUCCGAGCGGUACCCAAGUCCUUGUGCAACGUGAUAGUU"
    # s1 =       "..........(((((.(((....(((((((.(((.((.((.(((..((((((......(((....)))((((.((.....((((........)))))).))))....))))))))))).)).....))).)))))))....))))))))."
    # s2 =       "..........((((((((....(.((((((.((.....(((.((((((((.(((....)))))).....((((((((...))).))))).)))))..)))..(.((((((.(...))))))).)...)).)))))))...))).)))))."

    # Ex 139
    sequence = "AGUUCAGCUACUUACACUGCGGCCUACGGGGGAAUCACUGACCGCGUUAUACCCUAUAUUCGUAUAUGCACGUACAAUACCAACAGUAGGAAUACUUUCUGGGGAGUAGAGUAGAAAUCUUGGCCCCCGACUUGACUGUUAAUCCUCUAC"
    s1 =       ".((...(((...........)))..))(((((..(((......(((.(((((.........))))))))...(((..((((..(((.(((....))).)))..).)))..))).......))))))))(((......))).........."
    s2 =       "...........(((....(((((...((((((..(((......((.((((.(((((.((((.((..((...(((...)))...)).)).))))......))))).)))).))........)))))))))....).)))))))........"
    section_1 = [1, 36, 121, 142]
    section_2 = [36, 121]

    # merge ideas txt
    # sequence = "AACGCCGCGCCGGCCAGGCUGCCAGUUGUCUUAUUGAGUUGGCCAUCUACGUAAUACUUUUGUAGCCGACGACUAUAACG"
    # s1 =       "(((...(((((.....)).)))..)))...((((...((((((....((((.((....))))))))))))....)))).."
    # s2 =       "...((.(((((.....))).))..)).(((.......((((((....((((.........)))))))))))))......."
    # section_1 = [0,27]
    # section_2 = [27,78]

    # will always fail, not possible
    sequence = "GCCGUGUGCUCUGGCUCUAGACUCUAACACAGGGCGGUUGUCCCUAUUACUGGUACUCGU"
    s1       = "(((.((((....(........).....)))).)))....((.((.(....))).))...."
    s2       = "(((.((((..(((....))).......)))).)))((.((.((........))))))..."
    section_1 = [0,35]
    section_2 = [35,59]


    search_width=128
    
    new=True
    Debug = False
    Verbose = True


    paths_1 = pathfinder_cpp.find_path(sequence, s1, s2, section=section_1, search_width=search_width, Debug=Debug, Verbose=Verbose)
    paths_2 = pathfinder_cpp.find_path(sequence, s1, s2, section=section_2, search_width=search_width, Debug=Debug, Verbose=Verbose)



    Debug = False

    merge_search_width = 0.2

    start = time.time()
    merge_findpath(sequence, s1, s2, paths_1, paths_2, [],[], Debug=Debug, Verbose=Verbose, merge_search_width=merge_search_width)
    print ("time elapsed", time.time()-start)




    # sequence = "CCCUCAAUCGGUUCUUAAGACGUACUGNNNCCCGUUGGUAUGAAGAAAUUUGCUGGGAGAAAA"
    # s1 =       "..(((..(((((........(((((((...).....)))))).........))))))))...."
    # s2 =       ".((.(((....(((((....(((((((...).....)))))))))))..)))..))......."  
    # pathfinder_python.find_path(sequence, s1, s2, Verbose=True)    
    # pathfinder_python.find_path(sequence, s1, s2, Verbose=True)

    # sequence = "CCCUCAAUCGGUUCUUAAGACGUACUGCGCGUUUCACCAGACCANNNNNUGGUAUGGGAGGCGUGCCCGUUGGUAUGAAGAAAUUUGCUGGGAGAAAA"
    # s1 =       "..(((..(((((........(((((((......................................).....)))))).........))))))))...."
    # s2 =       ".((.(((....(((((....(((((((......................................).....)))))))))))..)))..))......."  
    # s1, s2 = s2, s1
    # r = pathfinder_python.find_path(sequence, s1, s2, Verbose=True)
    # print (r)


    # sequence = "CCCUCAAUCGGUUCUUAAGACGUACUGCGCGUUUCACCAGACCANNNNNUGGUAUGGGAGGCGUGCCCGUUGGUAUGAAGAAAUUUGCUGGGAGAAAA"
    # s1 =       "..(((..(((((........(((((((((((((((.(((.((((.....)))).)))))))))))).....)))))).........))))))))...."
    # s2 =       ".((.(((....(((((....(((((((((((((((((((((......))))))....))))))))).....)))))))))))..)))..))......."  

    # section = (22,75)
    # print ("orig. pathfinder: sE =", pathfinder.pathfinder(sequence, s1, s2, search_width=1000, section=section, verbose=True).sE)
    # print (pathfinder.pathfinder(sequence, s1, s2, search_width=1000, section=section, verbose=True).moves_pos)
    
    # section = (1, 22, 75, 94)
    # start = time.time()
    # Verbose = True
    print (pathfinder.pathfinder(sequence, s1, s2, search_width=search_width, verbose=Verbose).max_en)
    # # print (pathfinder.pathfinder(sequence, s1, s2, search_width=1000, verbose=True).moves_pos)
    # print ("time elapsed", time.time()-start)
    

    

