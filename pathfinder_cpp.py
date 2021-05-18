#!/usr/bin/env python3
# coding: utf-8

import time
import subprocess
import numpy as np

import pathfinder

import RNA
import helper

import networkx as nx
# from networkx.drawing.nx_agraph import graphviz_layout

import matplotlib.pyplot as plt

from helper import path_class


def levenshtein(a,b):
    # from stackoverflow etc...
    "Calculates the Levenshtein distance between a and b."
    n, m = len(a), len(b)
    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        a,b = b,a
        n,m = m,n
        
    current = range(n+1)
    for i in range(1,m+1):
        previous, current = current, [i]+[0]*n
        for j in range(1,n+1):
            add, delete = previous[j]+1, current[j-1]+1
            change = previous[j-1]
            if a[j-1] != b[i-1]:
                change = change + 1
            current[j] = min(add, delete, change)
            
    return current[n]

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


# def plot_graph(G, filename, x=20, y=12):

#     G = G.to_undirected()
#         # visualization below
#     # plt.figure(figsize=(12, 8), dpi=250, facecolor='w', edgecolor='k')
#     plt.figure(figsize=(x, y), dpi=250, facecolor='w', edgecolor='k')
#     pos = nx.spring_layout(G, dim=2, iterations=50)
#     edge_labels = nx.get_edge_attributes(G, "label")
#     node_labels = nx.get_node_attributes(G, "label")
#     # label offset
#     pos_attrs = {}
#     for node, coords in pos.items():
#         pos_attrs[node] = (coords[0], coords[1] - 0.04)

#     formatted_edge_labels = {(elem[0],elem[1]):edge_labels[elem] for elem in edge_labels} # use this to modify the tuple keyed dict if it has > 2 elements, else ignore

#     colors = [node[1]['color'] for node in G.nodes(data=True)]
#     nx.draw(G, pos, node_color=colors)
#     nx.draw_networkx_labels(G,pos_attrs,node_labels,font_size=14,font_color='b')
#     nx.draw_networkx_edge_labels(G,pos,edge_labels=formatted_edge_labels,font_color='red')

#     # plt.show()
#     plt.savefig(filename)




def find_path(sequence, s1, s2, search_width=500, section=None, Debug=False, Verbose=False, plot_graph=False):
    


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

    # findpath_cmd = "/mnt/c/Users/Max/Google\ Drive/arbeit/cppfinder/cppfinder/FINDpath"
    # findpath_cmd = "/mnt/c/Users/Max/Google\ Drive/arbeit/cppfinder/cppfinder/findpath_f"
    # findpath_cmd = "/mnt/c/Users/Max/Google\ Drive/arbeit/cppfinder/cppfinder/findpath_f_debug"
    # findpath_cmd = "/mnt/c/Users/Max/Google\ Drive/arbeit/cppfinder/cppfinder/findpath_f_2" # additional output
    # findpath_cmd = "/mnt/c/Users/Max/Google\ Drive/arbeit/cppfinder/cppfinder/findpath_f_min"
    # findpath_cmd = "/mnt/c/Users/Max/Google\ Drive/arbeit/cppfinder/cppfinder/findpath_f_minmax"
    findpath_cmd = "/mnt/c/Users/Max/Google\ Drive/arbeit/cppfinder/cppfinder/findpath_f_minmax2"
    # findpath_cmd = "/mnt/c/Users/Max/Google\ Drive/arbeit/cppfinder/cppfinder/findpath_f_minmax_debug"
    # findpath_cmd = "/mnt/c/Users/Max/Google\ Drive/arbeit/cppfinder/cppfinder/findpath_orig"
    # findpath_cmd = "/mnt/c/Users/Max/Google\ Drive/arbeit/cppfinder/cppfinder/findpath_orig_th"


    # 60.0 works with findpath_f_minmax but not with minmax2...

    # findpath_cmd = "/scratch/maxf/cppfinder/findpath_f"
    # findpath_cmd = "/scratch/maxf/cppfinder/findpath_f_debug"
    # findpath_cmd = "/scratch/maxf/cppfinder/findpath_f_all"



    if Debug:
        print ("input")
        print (sequence)
        print (s1)
        print (s2)


    start_findpath = time.time()

    cmd = f'printf "{sequence}\n{s1}\n{s2}" | {findpath_cmd} -m {search_width} -v'

    # cmd = f'echo -e "{sequence}\n{s1}\n{s2}" | {findpath_cmd} -m {search_width} -v'

    # print (cmd)

    result = subprocess.check_output(cmd, shell=True, encoding="utf8")

    end_findpath = time.time()-start_findpath

    if Verbose:
        print (result.split("---")[-1])
        # print (result)

    if Debug:
        print (result.split("---")[:-1])


    max_en_list = []

    # print (result)
# 
    sE = float(result.split()[-1])
    # # return without postprocessing
    # sE = float("-inf")
    # for line in result.split():
    #     if len(line)!=len(s1) and line[0] != "m":
    #         e = float(line)
    #         if e>sE:
    #             sE = e
    #     # if line[0] == "m":
    #     #     sE = float(line.split(',')[0].split(':')[1])/100
    # return_path = path_class()
    # return_path.max_en = sE
    # return return_path



    paths = result.split("---")[1].split("\n")
    paths = [path for path in paths if path.strip()] # remove newlines "" at start and end
     
    # 

    orig_path = paths[-1][4:].replace(' ', '').replace(')', '').replace('(', '').split(";")
    
    # paths = paths[:-1]
    # paths = list(set(paths))

    # print (paths)
    # print (return_path)

    # original output findpath
    fc = RNA.fold_compound(sequence)    
    bp_dist = RNA.bp_distance(s1, s2)
    current_ptable = list(RNA.ptable_from_string(s1))
    end_ptable = list(RNA.ptable_from_string(s2))
    

    return_path = path_class()

    return_path.e1_en = fc.eval_structure_pt(current_ptable)/100
    return_path.e2_en = fc.eval_structure_pt(end_ptable)/100
    return_path.dist = bp_dist    
    return_path.max_en = sE
    return_path.s1 = s1
    return_path.s2 = s2


    G = nx.DiGraph()


    current_ptable = list(RNA.ptable_from_string(s1))
    e = fc.eval_structure_pt(current_ptable)/100
    last_node = p_to_s(current_ptable)+str(0)
    first_node = last_node
    G.add_node(last_node)
    G.nodes[last_node]['color'] = "g"
    G.nodes[last_node]['label'] = str(0)+" / "+str(e)

    # paths.sort(key=len, reverse=True)
    counter_fwd = 0    
    counter_bwd = 0

    paths_info = []

    # process all paths
    for k, path in enumerate(paths):

        #workaround bwd
        temp_pt2 = list(RNA.ptable_from_string(s2))
        e_init = fc.eval_structure_pt(temp_pt2)/100

        mode = -1
        if (path[0:3]) == "FWD":
            # path_list = [(0,0,return_path.e1_en)]
            mode = True
        elif (path[0:3]) == "BWD":
            # path_list = [(0,0,return_path.e2_en)]
            mode = False
        else:
            continue

        path_list = []

        path = path[4:].replace(' ', '').replace(')', '').replace('(', '').split(";")
        
        max_en = float('-inf')

        # assumes we may get incomplete paths from pathfinder
        if len(path) > bp_dist: 
            if mode == True: # FWD path
                for l, move in enumerate(path):
                    if move == '': # end reached
                        continue

                    move = move.split('/')
                    it = int(move[0])
                    i = int(move[1])
                    j = int(move[2])
                    e = int(move[3])/100       
                    if i==0 and j==0: 
                        continue # end reached

                    path_list.append((i,j,e))

                
            if mode == False: # BWD

                for l, move in enumerate(reversed(path)):
                    if move == '': # end reached
                        continue
                    # print (move)
                    move = move.split('/')
                    it = int(move[0])
                    i = int(move[1])
                    j = int(move[2])
                    e = int(move[3])/100       

                    # path_list.insert(0,(i,j,e))
                    path_list.insert(0,(i,j,e_init))
                    e_init = e
                    



        
        path_list.insert(0,(0,0,return_path.e1_en)) # last step, dont interfere with bwd buildup
        
        # print ("current path:", path_list)

        opt_list = [x[2] for x in path_list] # energy list [kcal/mol, kcal/mol, ...]
        max_en = max(opt_list)
        max_en_list.append(max_en)
        max_en_pos = np.argmax(opt_list)

        # if max_en != -1.4:continue

        len_before = len(opt_list[0:max_en_pos])
        len_after = len(opt_list[max_en_pos+1:])

        if opt_list[0:max_en_pos] == []:
            num_before = 1
            min_before = opt_list[0]
            # min_before = [opt_list[0]]
            max_before = opt_list[0]
            opt_before = [opt_list[0]]
        else:
            num_before = len(opt_list[0:max_en_pos])
            min_before = np.min(opt_list[0:max_en_pos])
            # min_before = sorted(opt_list[0:max_en_pos])
            max_before = np.max(opt_list[0:max_en_pos])
            # opt_before = [sorted(opt_list[0:max_en_pos])[0]] + sorted(opt_list[0:max_en_pos], reverse=True)[:-1]
            opt_before = sorted(opt_list[0:max_en_pos], reverse=True)

        
        if opt_list[max_en_pos+1:] == []:
            num_after = 1
            min_after = opt_list[-1]
            # min_after = [opt_list[-1]]
            max_after = opt_list[-1]
            opt_after = [opt_list[-1]]
        else:
            num_after = len(opt_list[max_en_pos+1:])
            min_after = np.min(opt_list[max_en_pos+1:])
            # min_after = sorted(opt_list[max_en_pos+1:])
            max_after = np.max(opt_list[max_en_pos+1:])
            # opt_after  = [sorted(opt_list[max_en_pos+1:])[0]] + sorted(opt_list[max_en_pos+1:], reverse=True)[:-1]
            opt_after  = sorted(opt_list[max_en_pos+1:], reverse=True)


        # t_before = min_before/len_before
        # t_after  = min_after/len_after
        # t_before = len_before
        # t_after  = len_after
        t_before = mode
        t_after  = not mode

        num = num_before/num_after

        # opt_after = sorted(opt_before+opt_after, reverse=True)

        
        

        # print ("opt", opt_list, max_after)

        # sort by: we need the fwd & bwd path with lowest max_energy and secondly lowest min_before

        # print("opt", max_en, max_en_pos, min_before)
        return_path.paths.append(path_list)
        paths_info.append((max_en, min_before, max_before, min_after, max_after, t_before, t_after,opt_before, opt_after, num, opt_list, path_list))


    # sort all 
    # check similarity of results:


    # print (return_path.paths)

    def add_nodes(path, Verbose=False):

        current_ptable = list(RNA.ptable_from_string(s1))
        last_node = p_to_s(current_ptable)+str(0)
        e = fc.eval_structure_pt(current_ptable)/100

        # if path[-1][2][0] != -18:
        #     return

        plot_paths.append(path[-1])
        # print (path[0:-1])
        # if Verbose:
            # print (f"B: {0:2},{0:2} / {p_to_s(current_ptable)} {round(e,2):6.2f} ({0:4}, {0:4})")

        for l,(i,j,en) in enumerate(path[-1]):
            l = l-1 # this is excellent code
            if i==0 and j==0:
                continue
            if i<0:
                current_ptable[-i] = 0
                current_ptable[-j] = 0
            else:
                current_ptable[i] = j
                current_ptable[j] = i  

            e = fc.eval_structure_pt(current_ptable)/100

            current_node = p_to_s(current_ptable)+str(l+1)

            # print(current_node)

            G.add_node(current_node)
            G.nodes[current_node]['color'] = "r"
            G.nodes[current_node]['label'] = str(l+1)+" / "+str(e)
            G.add_edge(last_node, current_node)
            G[last_node][current_node]['label'] = (i,j)
            G[last_node][current_node]['weight'] = 1

            last_node = current_node

            old_s = p_to_s(current_ptable)
            # if Verbose:
            #     print (f"B: {k:2},{l+1:2} / {old_s} {round(e,2):6.2f} ({i:4}, {j:4})")

        if Verbose:
            # print ("max_en:", path[0], path[1], path[3], path[6], path[7], path[8], path[9])
            print ("max_en:", path[0], "moves:", path[-1])
      
            print ("~~~")
        # print ("saddle:")


    def find_next_path(sort_function):
            # check if this path is just a worse version of what we already have?
   # 

        while True:
            if paths_info==[]:
                break

            paths_info.sort(key=sort_function)
            current_path = paths_info[0]
            current_en_list = paths_info[0][-2]            

            exists = False
            for other_path in plot_paths:
                en_list = [x[2] for x in other_path] # energy list [kcal/mol, kcal/mol, ...]      

                # print (en_list, levenshtein(current_en_list,en_list))          
                if levenshtein(current_en_list,en_list)==0:
                    # path already exists
                    # print (en_list, levenshtein(current_en_list,en_list))
                    exists = True
            # break
            if exists:
                # print (current_en_list, "ignore")
                paths_info.pop(0)
                # try again
            if not exists:
                # print (current_en_list, "add")
                # add_nodes(paths_info[0], Verbose=True)
                add_nodes(paths_info[0], Verbose=False)
                paths_info.pop(0)
                break



    # add nodes to G

    plot_paths = []

    paths_info.sort(key=lambda x: (x[0], x[1], x[2]))
    return_path.paths = [(i[0],i[1]) for i in paths_info[0][-1]]
    # print ("test, ", paths_info[0][-1])

    find_next_path(lambda x: (x[0], x[1], x[2]))
    find_next_path(lambda x: (x[0], x[3], x[4]))

    find_next_path(lambda x: (x[6], x[1], x[2], x[0]))
    find_next_path(lambda x: (-x[6], x[3], x[4], x[0]))

    find_next_path(lambda x: (x[1], x[7]))
    find_next_path(lambda x: (x[3], x[8]))

    find_next_path(lambda x: (x[7]))
    find_next_path(lambda x: (x[8]))






    # print ("rest:", len(paths_info))    
    # if len(paths_info)==38:
    #     for path in paths_info:        
    #         add_nodes(path)
    #     add_nodes(paths_info[34])

    # print ("inner outer", len(paths_info))

    # if len(paths_info)==20:
    #     for path in paths_info:        
    #         add_nodes(path, Verbose=True)


    # for path in paths_info:        
    #     add_nodes(path)

    



    return_path.bp_dist = bp_dist

    return_path.max_en = min(max_en_list)

    return_path.runtime = end_findpath

    return_path.runtime_ext = time.time()-start_findpath-return_path.runtime


    # argmin here, find best path etc...


    current_ptable = list(RNA.ptable_from_string(s2))
    e = fc.eval_structure_pt(current_ptable)/100
    last_node = p_to_s(current_ptable)+str(len(orig_path)-1)


    # G.add_node(last_node)

    G.nodes[last_node]['color'] = "g"
    G.nodes[last_node]['label'] = str(bp_dist)+" / "+str(0)
    G.nodes[first_node]['color'] = "g"





    # print("done", bp_dist, counter_fwd, counter_bwd)

    return_path.G = G




    return return_path





if __name__ == '__main__':
    section = None

    # >sbox_leader
    # sequence = "AAAUGAUUCAACUUUUCAAAAAAUUAAUAACAUUUUCUCUUAUCGAGAGUUGGGCGAGGGAUUGGCCUUUUGACCCCAACAGCAACCGACCGUAAUACCAUUGUGAAAUGGGGCGCACUGCUUUUCGCGCCGAGACUGAUGUCUCAUAAGGCACGGUGCUAAUUCCAUCAGAUUGUGUCUGAGAGAUGAGAGAGGCAGUGUUUUACGUAGAAAAGCCUCUUUCUCUCAUGGGAAAGAGGCUUUUUGU"
    # s1 =       "............................((((((((((((((((....((((((((((((......))))))..)))))).(((.(((.((......(((((....)))))(((((..........)))))(((((....)))))....))..))))))........((((((...))))))..))))))))))..)))))).......((((((((((((((((....)))))))))))))))).."
    # s2 =       "............................(((((((((((((((((((.((((((.(((((.....)))))....))))))(((.((((.((......(((((....)))))(((((..........)))))(((((....)))))....))..)))))))..)))..((((((...))))))..)))))))))))..))))).......((((((((((((((((....)))))))))))))))).." 
  
  
    # sequence = "AAAUGAUUCAACUUUUCAAAAAAUUAAUUUUUCUCUUAUCGAGAGCGAGGGAUUGGCCUUUUGACCAGCAAUGUGAGUGCUAAUUCCAUUGUAGAGAUGAGAGAGGCAUUUACGUACUCAUUGU"
    # s1 =       "............................((((((((((((....(((((((......))))))..).(((.(....))))........(...)..))))))))))..)).......(....).."
    # s2 =       "............................(((((((((((((((.(.(((((.....)))))....)(((.((....)))))..)))..(...)..)))))))))))..).......(....).."

    # sequence = "ACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUAAGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUU"
    # s1 =       "..(((((((.....)))))))...(((((((((((((.......))))))).)))))).((((((((((.....))))))))))."
    # s2 =       "..(((((((.....)))))))................(...(((((....)))))...)((((((((((.....))))))))))."
    # sequence = "ACGAUUUCUAACGAAUUUUUUAAGUGCUUCUUGCUUAAGCAAGUUUCCCCUCAUU"
    # s1 =       "..(.....)...(((((((((((((.......))))))).)))))).(.....)."
    # s2 =       "..(.....)................(...(((((....)))))...)(.....)."


    # sequence = "AUGAACAGAUGGUACCUCGCGGCGGGGCCUGACCCAUCUACAGUUUGUGCAG"
    # s1 =       "((((((((((((..(((((...)))))......))))))...))))))...."
    # s2 =       ".((..((((..(((((....))..(((.....)))...)))..))))..))."


    # sequence = "AACGCCGCGCCGGCCAGGCUGCCAGUUGUCUUAUUGAGUUGGCCAUCUACGUAAUACUUUUGUAGCCGACGACUAUAACG"
    # s1 =       "(((...(..............)..)))...((((...((((((....((((.((....))))))))))))....)))).."
    # s2 =       "(((...(..............)..)))..........((((((....((((.........)))))))))).........."

    # sequence = "GGGCAGUAGUGUCCUACUGCCCAACCCGGUAGAGGGAUCGAGUGGUCUAGGCGAAUGAUACGAGAGGCUCACUGGGGACGAUGGGAGUGAUCUACUCGAUGUUGCGUCCGAUCACAUCCA"
    # s1 =       "(((((((((....)))))))))..((((((.(((...((..(((.((.........)))))..))..))))))))).....((((.((((((....((......))...)))))).))))"
    # s2 =       "(((((((((....)))))))))..(((......)))...((((((((..((((....((.(((((((.(((((.....(...)..)))))))).))))))....)))).)))))).)).."


    # sequence = "UGGAGGUAGAAGCGUUAGAGCAAUACUUCUACAGAGACCACGUGAGGUAGCUGUCUUAAGUCGUUGCUUACUGCGUGAGAUUGCAGCCAAUUAGCGGCUGUCUCAUCCGUACAUAAGGUA"
    # s1 =       ".(((((((...((......))..)))))))((((..(((......)))..)))).................((((((((((.((.((......)).)).)))))))..)))........."
    # s2 =       "..(((((((.(.(((.((((.....(((.(((........)))))).....))))..)))))))))))...........(((((.(.((((((..(...)..)))).))).)))))...."

    # # 60.1 inner
    sequence = "GCCAACAAACCGUGAUGGGCUAUGUUUGAUUCAUCCUAUUUAUGUUUUUCGAGAUGCGCG"
    s1       = "...........(((..(((..((.....))...)))....)))................."
    s2       = "...........((((((((..(((.......)))))))).)))................."

    # # 60.2 outer
    sequence = "UGAAGACCCAUUGAGUAACGACACCGCACGGCGCAUGGCGUCAGAGUAGCACUGCCUCGU"
    s1       = "....(((................................)))((.((......))))..."
    s2       = "....(((................................))).(((..((...))))).."
    # # 60.3 outer
    # s1       = ".......((((.(.(...((....)).....).)))))......................"
    # s2       = ".......((((...((..((........))..))))))......................"

    #hok/sok
    # sequence = "GGCGCUUGAGGCUUUCUGCCUCAUGACGUGAAGGUGGUUUGUUGCCGUGUUGUGUGGCAGAAAGAAGAUAGCCCCGUAGUAAGUUAAUUUUCAUUAACCACCACGAGGCAUCCCUAUGUCUAGUCCACAUCAGGAUAGCCUCUUACCGCGCUUUGCGCAAGGAGAAGAAGGCCAUGAAACUACCACGAAGUUCCCUUGUCUGGUGUGUGUUGAUCGUGUGUCUCACACUGUUGAUAUUCACUUAUCUGACACGAAAAUCGCUGUGCGAGAUUCGUUACAGAGACGGACACAGGGAGGUGGCGGCUUUCAUGGCUUACGAAUCCGGUAAGUAGCAACCUGGAGGCGGGCGCAGGCCCGCCUUUUCAGGACUGAUGCUGGUCUGACUACUGAAGCGCC"
    # s1       = "(((((((.((.(((((((((.(((((((....(((((....))))).))))))).)))))))))..(((.(((.(((.((..((((((....)))))).)).))).)))))).....(((.((.(((((((((....(((((((.((((((...))))..))..))).))))(((((((((.((((....(((((((((((((((.((((((((.(((((...))))).)))))))).)))...((((..(((((..((((...))))..)))))..))))..)))))).)))))))))).)).)))))))((.(((...........)))))..(((((((((((((....))))))))..))))).)))))).))).)))))..)).)))))))"
    # s2       = "((.(((.....(((((((((.(((((((....(((((....))))).))))))).))))))))).....))).)).......((((((....))))))......(((((((((..((((.......))))..)))).))))).....((((...))))..(((.....(((((((((((((.((((....(((((((((((((((.((((((((.(((((...))))).)))))))).)))...((((..(((((..((((...))))..)))))..))))..)))))).)))))))))).)).))))))))))).....))).((.....))......(((((((((....)))))))))..................................."


    # sequence = "UGCUAGUCAGGUUGUCAGUCUGUUUUUCGACCCGGACCCCGUAAAGACCCUGGCAAACAUUUCCCCGAAACACGUCGCGGGACGUAAGUCUGAGGCGCCAUGCGAAAACCAGACCGCUGGAGCGACGACCUCAGGAAACUAUGGAGACUU"
    # s1       = "((((((...((((..(.(((((..........)))))...)....))))))))))......(((((((......))).))))((((..(((((((((...(((.....((((....)))).))).)).)))))))....))))......."
    # s2       = "....((((.(..((((((...((((((.....(((...))).)))))).))))))..)...(((((((......))).))))((((.((((((((((((..(((..........)))..).)))....))))))...))))))..))))."
    #
    # sequence = "UGCUAGUCAGGUUGUCAGUCUGUUUUUCGACCCGGACCCCGUAAAGACCCUGGCAAACAUUAAACACGAAAGUGCGCCAUGCGAAAACCAGACCGCUGGAGCGACGAGGAAACGGAGACUU"
    # s1       = "((((((...((((..(.(((((..........)))))...)....))))))))))......(......)(..((((...(((.....((((....)))).))).)).))....)......."
    # s2       = "....((((.(..((((((...((((((.....(((...))).)))))).))))))..)...(......)(.(((((((..(((..........)))..).)))....)...)))..))))."

    # s1       = "...................................................................................(((..(((((((((...(((.....((((....)))).))).)).)))))))....)))........"
    # s2       = "...................................................................................(((.((((((((((((..(((..........)))..).)))....))))))...)))))........"


    # sequence = "AUGUUCUCAGUUGUUGUUGGAGCAGUCCCGUAGGAGGUCCCCAUUGGGCCCUUUUCUAGGGCUGCGAAGUCGGGCUCCAG"
    # s1       = ".((....)).......((((((((((((....((((((((.....))))).)))....))))))((....))..))))))"
    # s2       = "...........((..(((.((((((.(((.(((((((.(((....))).)))..)))))))))))....)).)))..))."
    
    # section = (1, 24, 62, 80)

    # sequence = "AUAGGUUGUGUGGAAGUCGACGGCUCUAGCGGAGGUUGCAAUCCGAGGAUGCUACGGAAUUAACGGGUUAAAGUGAAGCCGCGUCUAUCAAAGUUUUACCUUUUAGUUAGCCACGUCGCCGGUCGCGUAUAAGCCGGGUACCCGAGUACU"
    # s1       = "...(.(((((((((.(.(((((...(((((((((((.....((((((....)).))))((..(((((((.......)))).)))..)).........))))))..)))))...))))).)..)).))))))).)..(((((....)))))"
    # s2       = "...(((((((((((.(.((((((((((((...((((.(.(((..((((((((...((..((.((........)).)).)))))))).))...)))).)))).))))..))))..)))).)..)).))))))).)).(((((....)))))"
    # # section = [35, 98]
    # # section = [20, 115]
    # # section = [3, 20, 115, 150]

    # sequence = "AGGGGAAAGUUAUCGAACAGACCUUGUUUGGCUUAUGGCACUUGGGUAGAUUUUCUCUACCGUUGUCAUGAAAACGCACCGCAAGAGAAGUCACGAACCUGCGCCAGCGCACUACAGCUUGACUUUGGGAACUGAGUGGCGUAGUGACGG"
    # s1       = "..................................((((((....((((((.....))))))..))))))................................................................................."
    # s2       = "..................................((((((..((.(((((.....))))))).))))))................................................................................."
    # s1       = "...((...(((.(((((((.....))))))).((...................................)).)))...)).........(((((.....(((((((((........)))..((((.(....).))))))))))))))).."
    # s2       = "...((...((..(((.......(((((...((((...................................))....))...)))))........)))....)).))..((((((((.(((..((((........))))))))))))).))."
    # # section = [3, 34, 69, 149]

    # sequence = "CGACUUUCCUGAAUCCUGGCAGCUAACAUGGUGGCCGAGAGGGGCUCGCCCCCUCAACGAGUACCGAGAUUGCCUGCCUGCCGUGGCAAAUACGCGGGUCGCUCCGAGUACCCAACCUAGUAGUUAAUCGAUGCUCAAGCGAAAUUGUCG"
    # s1       = "..................((.((((...)))).))..................................................................................................................."
    # s2       = "..................(((.(......).).))..................................................................................................................."

    # sequence = "AAUAGAGCACUGUACCCAUCUCGUUACAUACUGUGUAACCGGAUUGUGUGUCCUGGUCCAAUAGAUAGGGUGCCGGUAGGCGGAACAUGGAAUGGGCUCGACACUUGUCACUCCGCGGCUCUGGUUAAGCCCCUCUCUGGGAGAUUGGGG"
    # s1       = ".................................................................................((((((.((..((.......)))))))...)))...................................."
    # s2       = ".................................................................................(((...(((((..(......)..)).))).)))...................................."
    
    # Ex50

    # sequence = "UCAUCCAAAGAACUCGUCAUGGGCGUCAACUUCAGUCGGUAAAGCUCACAGGUGACUUUGAUUUAAGCAAUGAGAUCUUGGACCUCUGCUAGAAGCUAUGACCUUGUCUGCCAAGCUUUUACUCAAAAAGACUCGUAUUCUCAAUCGUGA"
    # s2       = ".....(((.((.(((((.......(((((..(((.(..((.......)).).)))..))))).......))))).)))))......................................................................"
    # s1       = ".....(.((((.(((((.......(((.((((..((.((.....)).))))))))).............))))).)))))......................................................................"
    
    # Ex 115
    # sequence = "GAUAAGUAAAGGAAGCGCUUCCUCACGCGGACAGUGAUUCAGCAUAGCCAAACAACCUCUUAUUCAAGGUCUUUGGGUACUCGUCGUCUUGGAUCAAAUUAGUAGGCCAAAAAGGCGGUAGGUACUACGAUAACGAGAGGUUCUCUGCCA"
    # s1 =       "...........((..(((((((......))).))))..))....((.(((((...(((........)))..))))).))(((((..((.((((((..........(((.....)))....))).))).))..))))).(((.....)))."
    # s2 =       "...........((..(((((((....).)))..)))..)).(((.((.......(((((((.((.(..((.....((((((.((((((((((.((.........)))))...))))))).))))))))..))).))))))).)).))).."
    # section = [11, 17, 34, 149]
    # section = [17, 34]

    # sequence = "AACGCCGCGCCGGCCAGGCUGCCAGUUGUCUUAUUGAGUUGGCCAUCUACGUAAUACUUUUGUAGCCGACGACUAUAACG"
    # s1 =       "(((...(((((.....)).)))..)))...((((...((((((....((((.((....))))))))))))....)))).."
    # s2 =       "...((.(((((.....))).))..)).(((.......((((((....((((.........)))))))))))))......."
    # section = [0,27]
    # section = [27,78]

    # sequence = "AGGGGAAAGUUAUCGAACAGACCUUGUUUGGCUUAUGGCACUUGGGUAGAUUUUCUCUACCGUUGUCAUGAAAACGCACCGCAAGAGAAGUCACGAACCUGCGCCAGCGCACUACAGCUUGACUUUGGGAACUGAGUGGCGUAGUGACGG"
    # s1 =       "............................................................................................................................(.(....).)................"
    # s2 =       "............................................................................................................................(........)................"
    # section = [124, 134]



    # sequence = "UAAUUGUGUAUUAAUGAGCUUUAGGAGGGUUUUAUCACGUAUUUGUCUCUGAGCAAGAAG"
    # s1       = ".................((((..(((((....(((...)))....)))))))))......"
    # s2       = ".(((.....))).....((((.((((.((.(.........).)).)).))))))......"


    # sequence = "AGUUGUGCGGUAUAACGGGAGGUAGAAUGUAUUGGUGUAGACACUAUUGUCGAGCAUACA"
    # s2       = "....(((((......(....)......)))))...(((.(((......)))..)))...."
    # s1       = ".......(((((((((.....))....)))))))((((.((((....))))..))))..."


    # sequence = "GCUCUACCUGGAUGUAGUCGUGAAACGCGUGACUCCUCAGCUCGGAUGAAGGUACCAAUC"
    # s1       = "...........(.((...((((...))))..))).........................."
    # s2       = "...........(..((..(((...)))..))..).........................."

    # # lysC_B.subtilis
    # sequence = "GAAGAUAGAGGUGCGAACUUCAAGAGUAUGCCUUUGGAGAAAGAUGGAUUCUGUGAAAAAGGCUGAAAGGGGAGCGUCGCCGAAGCAAAUAAAACCCCAUCGGUAUUAUUUGCUGGCCGUGCAUUGAAUAAAUGUAAGGCUGUCAAGAAAUCAUUUUCUUGGAGGGCUAUCUCGUUGUUCAUAAUCAUUUAUGAUGAUUAAUUGAUAAGCAAUGAGAGUAUUCCUCUCAUUGC"
    # s1       = "..((((((..(.((((.((((........((((((...(((.......)))......)))))).......))))..)))))..(((((((((.(((.....))).)))))))))((((.((((((.....)))))).)))).((((((((....))))))))....))))))......((((((((((.....)))))))..)))...((((((((((.....))))))))))"
    # s2       = ".....(((..(.((((.((((........((((((...(((.......)))......)))))).......))))..)))))..(((((((((.(((.....))).)))))))))((((.((((((.....)))))).)))).((((((((....))))))))....)))(((((((((((((((((((.....)))))))..))...))))))))))................"


    # sequence,s1,s2 = helper.remove_constant_sections(sequence,s1,s2)
    print (sequence)

    # s2, s1 = s1, s2

    search_width = 500
    # Verbose = True
    Verbose = False
    # Debug = True
    Debug = False

    plot_graph=True
    # plot_graph=False

    start = time.time()
    r = find_path(sequence, s1, s2, section=section, search_width=search_width, Verbose=Verbose, Debug=Debug, plot_graph=plot_graph)
    
    print ("cppfinder:", r.max_en, round(r.runtime,2), round(r.runtime_ext,2), "s")
    # find_path(sequence, s1, s2, Verbose=False)

    # original implementation as comparison

    start = time.time()
    orig_Verbose = True
    orig = pathfinder.pathfinder(sequence, s1, s2, search_width=search_width, section=section, verbose=orig_Verbose).max_en
    print ("orig pathfinder:", orig, round (time.time()-start,4), "s")

    # find_path(sequence, s2, s1, Verbose=True)