#!/usr/bin/env python3
# coding: utf-8

import time
import subprocess
import sys

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

from dataclasses import dataclass

from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster import AgglomerativeClustering

import hdbscan

import findpath
import pathfinder
import pathfinder_i_cache


import RNA
import RNAxplorer


# @dataclass(frozen=True) # autohash
@dataclass
class structure_node:
    __slots__ = ('__dict__', 'structure', 'bp_dist_s1', 'bp_dist_s2', 'en')
    structure: str
    bp_dist_s1: int
    bp_dist_s2: int
    en: int

    max_en_s1: int = 0
    max_en_s2: int = 0

    def __str__(self):
        return f'{self.structure} -> {self.en:2.2f}'

    def __eq__(self, other):
        "Checks if equal."
        return self.structure == other.structure

    def __hash__(self):
        """Make the instances hashable."""
        return hash(self.structure)


# @dataclass(frozen=True) # autohash
# class intermediate:
#     structure: str
#     bp_dist_s1: int
#     bp_dist_s2: int
#     en: int
#     def __str__(self):
#         return f'{self.structure} -> {self.en:2.2f}'


def callRNAxplorer(seq, ref_struct1, ref_struct2, s1, s2, n=500):
    """
    Call RNAxplorer using ref_struct1, ref_struct2 as
    structures that are made equally probable within
    the ensemble. Tell RNAxplorer to generate n samples

    returns a Matrix2D object which contains all unique samples and the inputdata of the xplorer call.
    """
    RNA.cvar.uniq_ML = 1
    RNA.cvar.betaScale = 1

    # md = RNA.md("global")
    md = RNA.md()
    md.circ = 0
    md.uniq_ML = 1
    md.compute_bpp = 0
    md.betaScale = 1

    vc = RNA.fold_compound(seq, md, RNA.OPTION_MFE | RNA.OPTION_PF)
    extended_options = ""

    # samples = list (k,l,energy,[strucutres])
    tmpSamples = RNAxplorer.estimate_landscape(
        vc, ref_struct1, ref_struct2, n, extended_options)

    return_structures = set()

    # convert tmpSamples to a list with one structure per k,l-cell (min. energy for example).
    xplorerSamples = []
    for s in tmpSamples:
        uniqueStructures = list(set(s[3]))
        minEnergy = float("inf")
        minStructure = ""
        for structure in uniqueStructures:

            energy = int(RNA.energy_of_struct(seq, structure)*100)
            # return_structures.add((structure,energy))

            return_structures.add(structure_node(structure, RNA.bp_distance(structure, s1), RNA.bp_distance(structure, s2), energy))

            if energy < minEnergy:
                minStructure = structure
                minEnergy = energy
        xplorerSamples.append((s[0], s[1], minEnergy, minStructure))

    return list(return_structures)


def plot_graph(G, filename="g.png", x=12, y=8):

    G = G.to_undirected()

    plt.figure(figsize=(x, y), dpi=250, facecolor='w', edgecolor='k')
    pos = nx.spring_layout(G, dim=2, iterations=50)
    edge_labels = nx.get_edge_attributes(G, "weight")
    node_labels = nx.get_node_attributes(G, "label")

    # label offset
    pos_attrs = {}
    for node, coords in pos.items():
        pos_attrs[node] = (coords[0], coords[1] - 0.07)

    formatted_edge_labels = {(elem[0], elem[1]): edge_labels[elem] for elem in edge_labels}  # use this to modify the tuple keyed dict if it has > 2 elements, else ignore

    colors = [node[1]['color'] for node in G.nodes(data=True)]
    nx.draw(G, pos, node_color=colors)
    nx.draw_networkx_labels(G, pos_attrs, node_labels, font_size=15, font_color='b')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=formatted_edge_labels, font_size=11, font_color='red')

    # plt.show()
    plt.savefig(filename)


def add_node(G, i):
    
    s, dist1, dist2, eval = i.structure, i.bp_dist_s1, i.bp_dist_s2, i.en
    
    s_node = structure_node(s, dist1, dist2, eval)
    G.add_node(s_node)
    G.nodes[s_node]['color'] = "r"
    G.nodes[s_node]['label'] = (s_node.bp_dist_s1, s_node.bp_dist_s2)
    return s_node


def add_edge(G, n1, n2, weight):
    G.add_edge(n1, n2)
    G[n1][n2]['weight'] = weight



def detour(s, s1, s2):

    # fix: should also include missing bps 

    pt = RNA.ptable_from_string(s)
    pt1 = RNA.ptable_from_string(s1)
    pt2 = RNA.ptable_from_string(s2)


    non_direct_bp=0
    detour_i = []

    for i in range(1,pt[0]):

        if pt[i]>0 and i<pt[i]:
            if pt[i]!=pt1[i] and pt[i]!=pt2[i]:
                detour_i.append(i)
                # non_direct_bp += 1


    # print (s)
    # print (s1)
    # print (s2)
    # print (non_direct_bp)

    return detour_i
    # return non_direct_bp*1.0


def stabilize(s, detour_i, fc):

    s_eval = fc.eval_structure(s)

    pt = list(RNA.ptable_from_string(s))

    for i in detour_i:
        j = pt[i]
        pt[i] = 0
        pt[j] = 0
    
    s_new = RNA.db_from_ptable(pt)
    s_eval_new = fc.eval_structure(s_new)

    # print (s, s_eval)
    # print (s_new, s_eval_new)

    return s_eval_new-s_eval


def init(sequence, s1, s2):

    fc = RNA.fold_compound(sequence)
    s1_eval = fc.eval_structure(s1)
    s2_eval = fc.eval_structure(s2)

    search_width_multiplier = 2

    max_en = findpath.init_single_findpath(
        sequence, s1, s2, search_width_multiplier, True)

    bp_dist_total = RNA.bp_distance(s1, s2)

    print(max_en)

    G = nx.Graph()

    s1_node = add_node(G, structure_node(s1, 0, bp_dist_total, s1_eval))
    s2_node = add_node(G, structure_node(s2, bp_dist_total, 0, s2_eval))
    add_edge(G, s1_node, s2_node, max_en)

    # s3_node = add_node(G, s1, 0, bp_dist_total, s1_eval)
    # add_edge(G, s1_node, s3_node, max_en)

    # intermediates = callRNAxplorer(sequence, s1, s2, s1, s2, 4)
    # intermediates += callRNAxplorer(sequence, s1, s2, s1, s2, 2)
    # intermediates += callRNAxplorer(sequence, s1, s2, s1, s2, 2)

    # intermediates = [i[0] for i in RNAxplorer.repellant_sampling(fc)]
    # print ("here", RNAxplorer.repellant_sampling(fc))


    cmd = f'RNAxplorer -M RSH -n 200 --sequence {sequence}'
    result = subprocess.check_output(cmd, shell=True, encoding="utf8")
    result = result.split("\n")
    intermediates = []
    for line in result:
        line = line.split()
        if len(line) != 4: continue
        structure = line[1]
        energy = int(RNA.energy_of_struct(sequence, structure)*100)
        intermediates.append(structure_node(structure, RNA.bp_distance(structure, s1), RNA.bp_distance(structure, s2), energy))




    print ("RNAxplorer done")

    intermediates.sort(key=lambda x: x.en)


    x = [i.structure for i in intermediates]

    distance_s1s2 = np.zeros([len(x), 4])

    distances = np.zeros([len(x), len(x)])
    for i, ts1 in enumerate(x):
        for j, ts2 in enumerate(x):
            if ts1 == ts2:
                continue
            dist = RNA.bp_distance(ts1, ts2)
            distances[i][j] = dist
            distances[j][i] = dist
        
        distance_s1s2[i][0] = RNA.bp_distance(ts1, s1)
        distance_s1s2[i][1] = RNA.bp_distance(ts1, s2)
        distance_s1s2[i][2] = len(detour(ts1, s1, s2))
        # distance_s1s2[i][3] = int(fc.eval_structure(ts1)*100)

    print ([i.structure for i in intermediates])

    from sklearn.preprocessing import normalize
    from sklearn.preprocessing import MinMaxScaler

    # min_max_scaler = MinMaxScaler()
    # distance_s1s2 = min_max_scaler.fit_transform(distance_s1s2)    

    # print (distance_s1s2)

    # model = AgglomerativeClustering(
    # distance_threshold=0, n_clusters=None, affinity="precomputed", linkage="average", compute_full_tree=True)
    # model = model.fit(distances)

    # clusterer = AgglomerativeClustering(
    # distance_threshold=0, n_clusters=None, linkage="average", compute_full_tree=True)

    clusterer = AgglomerativeClustering(
    distance_threshold=5, n_clusters=None, linkage="average", compute_full_tree=True)
    clusterer.fit(distance_s1s2)

    # clusterer = hdbscan.HDBSCAN(metric='precomputed')    
    # clusterer.fit(distances)

    # clusterer = hdbscan.HDBSCAN()    
    # clusterer.fit(distance_s1s2)

    
    number_of_clusters = clusterer.labels_.max()

    print ("number of clusters:", number_of_clusters+1)


    for current_label in range(number_of_clusters+1):
        
        members = (clusterer.labels_==current_label).nonzero()[0]

        # print (current_label, "m:", members)

        min_en_cluster = float("inf")

        collect_s = []

        for i in members:


            s = intermediates[i].structure
            current_detour = detour(s, s1, s2)

            if len(current_detour)==0: continue

            collect_s.append(s)

            en = int(fc.eval_structure(s)*100)
            if en<min_en_cluster:
                min_en_cluster = en

            max_en_s1 = findpath.init_single_findpath(sequence, s, s1, search_width_multiplier, True)
            max_en_s2 = findpath.init_single_findpath(sequence, s, s2, search_width_multiplier, True)

            stabilize_en = stabilize(s, current_detour, fc)

            print (current_label, i, s, distance_s1s2[i][0], distance_s1s2[i][1], current_detour, "/", max(max_en_s1, max_en_s2), en, stabilize_en)

        print ("done1", current_label, "min_en", min_en_cluster)

        search_width = 128
        Debug = False
        Verbose = False

        # pathfinder_i_cache.find_path(sequence, [s1], collect_s, search_width=search_width, Debug=Debug, Verbose=Verbose)

        print ("done2", current_label, "min_en", min_en_cluster)

        # pathfinder_i_cache.find_path(sequence, [s2], collect_s, search_width=search_width, Debug=Debug, Verbose=Verbose)
        # break


    # intermediates = intermediates[0:50] # does not work

    for it, i in enumerate(intermediates):

        add_node(G, i)
        max_en_s1 = findpath.init_single_findpath(sequence, i.structure, s1, search_width_multiplier, True)
        add_edge(G, i, s1_node, max_en_s1)
        max_en_s2 = findpath.init_single_findpath(sequence, i.structure, s2, search_width_multiplier, True)
        add_edge(G, i, s2_node, max_en_s2)

        pathfinder.pathfinder(sequence, s1, i.structure, search_width=2*bp_dist_total, verbose=False)
        pathfinder.pathfinder(sequence, i.structure, s2, search_width=2*bp_dist_total, verbose=False)

        i.max_en_s1 = max_en_s1
        i.max_en_s2 = max_en_s2

        # print (i, max(i.max_en_s1, i.max_en_s2), i.max_en_s1, i.max_en_s2, "bp1", i.bp_dist_s1, "bp2", i.bp_dist_s2, clusterer.labels_[it])

    intermediates.sort(key=lambda x: (max(x.max_en_s1, x.max_en_s2), (x.bp_dist_s1+x.bp_dist_s2))  )


    print ("~~~")

    intermediates = intermediates[0:60]

    i_list = []
    # G.add_node(current_node)
    # G.nodes[current_node]['label'] = str(l+1)+" / "+str(e)
    for i in intermediates:
        print(i.structure, i.en, max(i.max_en_s1, i.max_en_s2), i.max_en_s1, i.max_en_s2, i.bp_dist_s1, i.bp_dist_s2)

        i_list.append(i.structure)

        # pathfinder.pathfinder(sequence, s1, i.structure, search_width=2*bp_dist_total, verbose=True)
        # pathfinder.pathfinder(sequence, i.structure, s2, search_width=2*bp_dist_total, verbose=True)
        i.max_en_s2 = findpath.init_single_findpath(sequence, i.structure, s2, search_width_multiplier, True)

        # break


    print (f"sequence = '{sequence}'")
    print (f"s1       = '{s1}'")
    print (f"s2       = '{s2}'")

    print (f"intermediates = {i_list}")


    # plot_graph(G)

    # G[last_node][current_node]['label'] = (i,j)
    # G[last_node][current_node]['weight'] = 1


if __name__ == "__main__":

    # sequence = "CACAUGGGAAGCAAUGUAGCUCUCUGAUUCAGCUCUGCUACCAAAUAAUUGGCCAGGUCCAGCGCGAUCUCUCGGUUGGGCAACUAUGUCAGUGGAGGGU"
    # s1 = "...((.((((((......))).))).))...........(((.....((((((...(((((((.(((....))))))))))......))))))....)))"
    # s2  = ".....(((.(((......)))))).......(((((.((((........((((..(((...((.(((((....))))).)).)))..)))))))))))))"

    sequence = 'UAAAAUGAUCACGGUUUCAGCUUUGGACAGGGCGUUCCACUAAACUCCUGGUGACAUAGAUAUAUUGGAUUGCAACUACUCGUCGGUCCGGUUGGCGUUC'
    s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'
    s2       = '....(((.(((((((....))).....(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........'

    # s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..' #s1
    # s2       = '....(((.((((...((((....))))(((((.(((......)))))))))))))))......((((((((((........).)))))))))........' #mfe





# ....(((.((((...((((....))))(((((.(((......)))))))))))))))......((((((((((........).)))))))))........ -> -1850.00 -540 -540 -1520
# ....(((.(((((((....))).....(((((.(((......)))))))))))))))......((((((((((........).)))))))))........ -> -1779.00 -540 -540 -1560

    # s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'
    # s2       = '..........(((...((((((..((((((((.(((......)))))))...(((...((....(((.....)))....))))).)))))))))))))..'

    # s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'
    # s2       = '..........(((...((((((..((((.(((.(((......))).)))...(((...((....(((.....)))....))))).)))))))))))))..'


    # simple 60 nt ex
    sequence = "GAAAGUUUUCAAGAAAGUAAUAUAGUGUCCCAGUCAGGGGGUAGUGGCGCGUUCUGGACGCAUAGUAUUC"
    s1 = "(((....))).........(((((((((.((((...(.(.((....)).).).)))))))).)).))).."
    s2 = ".................((.(((.((((.((((..(.(..((....)).).).))))))))))).))..."

    # 80 nt, S: -16.60 simple detour, -17.70 better one
    sequence = "CGGGAGCGAGUAUGGACCGGGUUGAAUUGGUACCUCGGCCGCCCUCUGGGGGUAGCCGACGGGCGGCUUCAGCCGGGCCC"
    s1 = ".............((.((.(((((((.(.((.(((((((.(((((...))))).))))).)))).).))))))).)).))"
    s2 = "((....)).....((.((.(((((((...........(((((((((.(.(.....)))).)))))))))))))).))))."

    # 80 nt, -7.20 direct, -8.30 simple detour, -8.40
    sequence = "CUGGGUCGUGGUGGCCUUUUAGAUACGAUUCACGAACGUAGCACGUUUCGGUCUCCGGAGACGCAAUGAUCUCGAGGGUA"
    s1 = ".(((((((((.(((....)))..))))))))).((.(((.((..(((((((...))))).)))).))).))........."
    s2 = ".(((((((((.............)))))))))(((((((...))))).)).((((.(((..........))).))))..."
    
    # 150 nt, -20.10 kcal/mol direct, indirect -24.50
    sequence = 'AAAAUAAUGUACCGGACAUUCGCGCACGACCACCAUAUGGCAGAGCAUGUGUCUGUGGACCCACUAUAGCUGGGGCGCUUAACCCCAGAAAAGUAUCUUCGGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAA'
    s1       = '.............((((((..((...(..(((.....)))..).))..))))))(((((((....(((.((((((.......)))))).....)))....)))))))(((......((.((.........)).))......)))......'
    s2       = '......(((.((((((...(((....))).(((((((((......))))))...)))((...(((....((((((.......))))))...))).)))))))).)))(((......((.((.........)).))......)))......'


    init(sequence, s1, s2)

    max_en_s1 = findpath.init_single_findpath(sequence, s1, s2, 2, True)
    print ("direct:", max_en_s1)

    # CACAUGGGAAGCAAUGUAGCUCUCUGAUUCAGCUCUGCUACCAAAUAAUUGGCCAGGUCCAGCGCGAUCUCUCGGUUGGGCAACUAUGUCAGUGGAGGGU
    # ...((.((((((......))).))).))...........(((.....((((((...(((((((.(((....))))))))))......))))))....))) [   0,    0 ] -24.60
    # ....(.((((((......))).))).)............(((.....((((((...(((((((.(((....))))))))))......))))))....))) [  -4,  -28 ] -23.50
    # ......((((((......))).)))..............(((.....((((((...(((((((.(((....))))))))))......))))))....))) [  -5,  -27 ] -23.00
    # .......(((((......))).))...............(((.....((((((...(((((((.(((....))))))))))......))))))....))) [  -7,  -25 ] -21.60
    # .......(.(((......)))..)...............(((.....((((((...(((((((.(((....))))))))))......))))))....))) [  -9,  -23 ] -19.10
    # .........(((......)))..................(((.....((((((...(((((((.(((....))))))))))......))))))....))) [  -8,  -24 ] -21.30
    # .......(.(((......)))).................(((.....((((((...(((((((.(((....))))))))))......))))))....))) [   8,   22 ] -20.70
    # ......((.(((......)))))................(((.....((((((...(((((((.(((....))))))))))......))))))....))) [   7,   23 ] -22.00
    # .....(((.(((......))))))...............(((.....((((((...(((((((.(((....))))))))))......))))))....))) [   6,   24 ] -24.50
    # ..(..(((.(((......)))))).).............(((.....((((((...(((((((.(((....))))))))))......))))))....))) [   3,   26 ] -22.90
    # ..((.(((.(((......)))))))).............(((.....((((((...(((((((.(((....))))))))))......))))))....))) [   4,   25 ] -24.80
    # ..((.(((.(((......))))))))..............((.....((((((...(((((((.(((....))))))))))......))))))....)). [ -40, -100 ] -24.00
    # ..((.(((.(((......))))))))...............(.....((((((...(((((((.(((....))))))))))......))))))....).. [ -41,  -99 ] -21.10
    # ..((.(((.(((......)))))))).....................((((((...(((((((.(((....))))))))))......))))))....... [ -42,  -98 ] -22.10
    # ..((.(((.(((......))))))))........(............((((((...(((((((.(((....))))))))))......))))))...)... [  35,   97 ] -17.40
    # ..((.(((.(((......))))))))........((...........((((((...(((((((.(((....))))))))))......))))))..))... [  36,   96 ] -19.90
    # ..((.(((.(((......))))))))........((.(.........((((((...(((((((.(((....))))))))))......)))))).)))... [  38,   95 ] -17.50
    # ..((.(((.(((......))))))))........((.((........((((((...(((((((.(((....))))))))))......))))))))))... [  39,   94 ] -20.10
    # ..((.(((.(((......))))))))........((.((.........(((((...(((((((.(((....))))))))))......))))).))))... [ -48,  -93 ] -17.50
    # ..((.(((.(((......))))))))........((.(((........(((((...(((((((.(((....))))))))))......))))))))))... [  40,   93 ] -19.70
    # ..((.(((.(((......))))))))........((.(((.........((((...(((((((.(((....))))))))))......)))).)))))... [ -49,  -92 ] -17.90
    # ..((.(((.(((......))))))))........((.((((........((((...(((((((.(((....))))))))))......))))))))))... [  41,   92 ] -21.80
    # ..((.(((.(((......)))))))).......(((.((((........((((...(((((((.(((....))))))))))......))))))))))).. [  34,   98 ] -22.30
    # ..((.(((.(((......))))))))......((((.((((........((((...(((((((.(((....))))))))))......)))))))))))). [  33,   99 ] -25.00
    # ..((.(((.(((......))))))))......((((.((((........((((....((((((.(((....))))))))).......)))))))))))). [ -57,  -81 ] -20.80
    # ..((.(((.(((......))))))))......((((.((((........((((...(((((((.(((....)))))))))...)...)))))))))))). [  57,   84 ] -20.20
    # ..((.(((.(((......))))))))......((((.((((........((((..((((((((.(((....)))))))))...))..)))))))))))). [  56,   85 ] -22.60
    # ..((.(((.(((......))))))))......((((.((((........((((..((.(((((.(((....))))))))....))..)))))))))))). [ -58,  -80 ] -21.00
    # ..((.(((.(((......))))))))......((((.((((........((((..((((((((.(((....))))))))...)))..)))))))))))). [  58,   83 ] -23.30
    # ..((.(((.(((......)))))))).....(((((.((((........((((..((((((((.(((....))))))))...)))..))))))))))))) [  32,  100 ] -24.50
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((.((((.(((....)))))))....)))..))))))))))))) [ -59,  -79 ] -20.40
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((..(((.(((....)))))).....)))..))))))))))))) [ -60,  -78 ] -17.50
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((...((.(((....)))))......)))..))))))))))))) [ -61,  -77 ] -16.60
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((....(.(((....)))).......)))..))))))))))))) [ -62,  -76 ] -14.60
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((......(((....)))........)))..))))))))))))) [ -63,  -75 ] -16.40
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((...(..(((....)))......).)))..))))))))))))) [  62,   81 ] -13.50
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((...((.(((....))).....)).)))..))))))))))))) [  63,   80 ] -16.20
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((...((.((......)).....)).)))..))))))))))))) [ -67,  -72 ] -15.20
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((...((..(......)......)).)))..))))))))))))) [ -65,  -74 ] -13.50
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((...((................)).)))..))))))))))))) [ -66,  -73 ] -16.50
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((...((.(............).)).)))..))))))))))))) [  65,   78 ] -18.10
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((...((.((..........)).)).)))..))))))))))))) [  66,   77 ] -18.70
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((...((.(((........))).)).)))..))))))))))))) [  67,   76 ] -20.60
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((...((.((((......)))).)).)))..))))))))))))) [  68,   75 ] -22.20
    # ..((.(((.(((......)))))))).....(((((.((((........((((..(((...((.(((((....))))).)).)))..))))))))))))) [  69,   74 ] -24.70
    # ..(..(((.(((......)))))).).....(((((.((((........((((..(((...((.(((((....))))).)).)))..))))))))))))) [  -4,  -25 ] -22.80
    # .....(((.(((......)))))).......(((((.((((........((((..(((...((.(((((....))))).)).)))..))))))))))))) [  -3,  -26 ] -24.40
    # S: -13.50 kcal/mol | B:  11.10 kcal/mol | E[start]:-24.60 E[end]:-24.40
