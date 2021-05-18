#!/usr/bin/env python3
# coding: utf-8

import sys
import random
import string

import numpy as np
import os
import subprocess
from dataclasses import dataclass, field

import matplotlib.pyplot as plt
from PIL import Image
from IPython.display import SVG, display

import RNA
import RNAxplorer
import pathfinder

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering
from sklearn.datasets import load_iris


def callRNAxplorer(seq, ref_struct1, ref_struct2, n=500):
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
        minEnergy = sys.float_info.max
        minStructure = ""
        for structure in uniqueStructures:

            return_structures.add(structure)
            energy = RNA.energy_of_struct(seq, structure)
            if energy < minEnergy:
                minStructure = structure
                minEnergy = energy
        xplorerSamples.append((s[0], s[1], minEnergy, minStructure))

    return list(return_structures)
    # return xplorerSamples

    # x = Matrix2D(seq, xplorerSamples, ref_struct1, ref_struct2)
    # return x

sequence = 'UAAAAUGAUCACGGUUUCAGCUUUGGACAGGGCGUUCCACUAAACUCCUGGUGACAUAGAUAUAUUGGAUUGCAACUACUCGUCGGUCCGGUUGGCGUUC'
s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'
s2       = '....(((.(((((((....))).....(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........'


s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..' #s1
s2       = '....(((.((((...((((....))))(((((.(((......)))))))))))))))......((((((((((........).)))))))))........' #mfe

s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..' #s1
s2       = '..........(((...((((((..(((((((..(((......))).)))..((((...((....(((.....)))....)))))))))))))))))))..' #A

s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..' #s1
s2       = '..........(((...((((((..((((..((.(((......))).))....(((...((....(((.....)))....))))).)))))))))))))..' #B

s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..' #s1
s2       = '..........(((...((((((..((((.(((.(((......))).)))(((..((..............))..)))........)))))))))))))..' #C
s1       = '..........(((...((((((..((((.(((.(((......))).)))(((..((..............))..)))........)))))))))))))..' #C
s2       = '..........(((...((((((..((((..((.(((......))).))....(((...((....(((.....)))....))))).)))))))))))))..' #B
s1       = '..........(((...((((((..((((..((.(((......))).))....(((...((....(((.....)))....))))).)))))))))))))..' #B
s2       = '....(((.((((...((((....))))(((((.(((......)))))))))))))))......((((((((((........).)))))))))........' #mfe



# sequence = 'CUGGGUCGUGGUGGCCUUUUAGAUACGAUUCACGAACGUAGCACGUUUCGGUCUCCGGAGACGCAAUGAUCUCGAGGGUA'
# s1       = '.(((((((((.............)))))))))(((((((...))))).)).((((.(((..........))).))))...'
# s2       = '.(((((((((.............)))))))))(((.(((.((...((((((...))))))..)).)))...)))......'

# sequence = 'GGAAGCCGGCGAGGCAGUACCAUUAUAUAGUUUGUCUUCCAAGAAUGGGUACGACCGCGGGACCGUUCGGUUAUCGUCUG'
# s1 = '((.((((((((.((..((((((((....((.....))......))))).)))..)).)).......)))))).)).....'
# s2 = '((.((((((.(((((((.((.........)))))))))))..((((((...((....))...)))))))))).)).....'

# s0 = '((.((((((((.((..((((((((....((.....))......))))).)))..)).)).......)))))).)).....'
# s1 = '.((((((((((.((..((((((((....((.....))......)))).))))..)).......))).))))).)).....'
# s1 = '.......(((((((..((((((((....((.....))......)))).))))..))....((((....)))).)))))..'

# s2 = '((.((((((.(((((((.((.........)))))))))))..((((((...((....))...)))))))))).)).....'

# s1 = '((..(((.....))).(((((((..(...............)..))).))))..))((((((((....)))).))))...' # mid
# s2 = '((.((((((.(((((((.((.........)))))))))))..((((((...((....))...)))))))))).)).....' # s2
# s2 = '((.((((((((.((..((((((((....((.....))......))))).)))..)).)).......)))))).)).....' # s1

# s1 = '((.((((((((.((..((((((((....((.....))......))))).)))..)).)).......)))))).)).....'
# s2 = '((..(((.....))).((((((((....((.....))......)))).))))..))(((.((((....))))..)))...'

# ....(.(((((.((..((((((((....((.....))......))))).)))..)).))...)))).............. 53 12 22 8 -11.2 -7.9
# (((((((.....))).(((....)))..........))))......((......))(((.((((....))))..)))... 2 21 44 18 -12.5 -5.2

# bp_dist0 = RNA.bp_distance(s1, s0)
bp_dist_all = RNA.bp_distance(s1, s2)

# s1 = '.......(((((((..((((((((....((.....))......))))).)))..))....((((....)))).)))))..'
# s2 = '((.((((((.(((((((.((.........)))))))))))..((((((...((....))...)))))))))).)).....'

fc = RNA.fold_compound(sequence)
sw = bp_dist_all*2
# sw = int(bp_dist_all/2)
x = callRNAxplorer(sequence, s1, s2, 2)

# for i in x:
#     print(i)


# distances = np.zeros([len(x), len(x)])

# for i, ts1 in enumerate(x):
#     for j, ts2 in enumerate(x):
#         if ts1 == ts2:
#             continue
#         dist = RNA.bp_distance(ts1, ts2)
#         distances[i][j] = dist
#         distances[j][i] = dist

# # np.set_printoptions(threshold=sys.maxsize)
# np.set_printoptions(threshold=1000)
# np.set_printoptions(linewidth=200)


# # find closest structure to s1 and s2
# bp_dist1 = 9999
# s1_closest = False
# bp_dist2 = 9999
# s2_closest = False

# for i, ts in enumerate(x):
#     dist1 = RNA.bp_distance(s1, ts)
#     if dist1 < bp_dist1:
#         bp_dist1 = dist1
#         s1_closest = (i, ts)
#     dist2 = RNA.bp_distance(s2, ts)
#     if dist2 < bp_dist2:
#         bp_dist2 = dist2
#         s2_closest = (i, ts)


# print(distances, sys.maxsize)
# print(s1_closest, bp_dist1)
# print(s2_closest, bp_dist2)

# # RNAxplorer.levelSaddlePoint(sequence, s1, s2, 2, 500, 100, 100)


# def plot_dendrogram(model, **kwargs):
#     # Create linkage matrix and then plot the dendrogram

#     # create the counts of samples under each node
#     counts = np.zeros(model.children_.shape[0])
#     n_samples = len(model.labels_)
#     for i, merge in enumerate(model.children_):
#         current_count = 0
#         for child_idx in merge:
#             if child_idx < n_samples:
#                 current_count += 1  # leaf node
#             else:
#                 current_count += counts[child_idx - n_samples]
#         counts[i] = current_count

#     linkage_matrix = np.column_stack([model.children_, model.distances_,
#                                       counts]).astype(float)

#     # Plot the corresponding dendrogram
#     dendrogram(linkage_matrix, **kwargs)
#     return

# def detect_local_minimum(fc, structure):
#     # perform gradient walk from sample to determine direct local minimum
#     pt = RNA.IntVector(RNA.ptable(structure))
#     fc.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
#     return RNA.db_from_ptable(list(pt))



# iris = load_iris()
# # X = iris.data

# # setting distance_threshold=0 ensures we compute the full tree.
# model = AgglomerativeClustering(
#     distance_threshold=0, n_clusters=None, affinity="precomputed", linkage="average", compute_full_tree=True)

# model = model.fit(distances)
# plt.title('Hierarchical Clustering Dendrogram')
# # plot the top three levels of the dendrogram
# plot_dendrogram(model, truncate_mode='level', p=9999)
# plt.xlabel("Number of points in node (or index of point if no parenthesis).")
# plt.savefig("dendogram.png", dpi=200)

# mem = {}


# def find_descendent(node):
#     N = len(x)
#     global mem
#     if node in mem:
#         return mem[node]
#     if node < N:
#         return [node]
#     pair = model.children_[node-N]
#     mem[node] = find_descendent(pair[0])+find_descendent(pair[1])
#     return mem[node]


# def find_ancestor(target):
#     N = len(x)
#     for ind, pair in enumerate(model.children_):
#         if target in pair:
#             return [target]+find_ancestor(N+ind)
#     return [ind+N]


# # ans = []
# # for a in find_ancestor(s1_closest[0])[::-1]:
# #     ans.append(find_descendent(a))

# # for a in ans:
# #     c = a[-1]
# #     print (c)

# x_ext = {i: x[i] for i in range(len(x))}

# # print (model.distances_)


# def join_s(s1, s2):

#     s3 = ""
#     for i, j in zip(s1, s2):
#         if (i == j):
#             s3 += i
#         elif i == ".":
#             s3 += j
#         elif j == ".":
#             s3 += i
#         else:
#             s3 += "x"
#     return s3


# for i, merge in enumerate(model.children_):
#     # print (i+len(x), merge, len(x))

#     # print (x_ext[merge[0]])
#     # print (x_ext[merge[1]])
#     # print (x_ext)

#     s3 = join_s(x_ext[merge[0]], x_ext[merge[1]])

#     x_ext[i+len(x)] = s3
#     # print (x_ext[i+len(x)], i+len(x))

#     # break

# # print (x_ext)

# # print (ans)

# # init
# init_se = pathfinder.pathfinder(sequence, s1, s2, search_width=sw).max_en
# print("init findpath:", init_se)

# last_descendent = set([s1_closest[0]])
# processed = set([s1_closest[0]])

# for i in find_ancestor(s1_closest[0]):

#     if i == s1_closest[0]:
#         continue
#     current_descendent = set(find_descendent(i))

#     # print (i, x_ext[i])
# #     mem = find_descendent(i)
# #     # print (mem)
#     # last_descendent = current_descendent-last_descendent

#     candidates = []

#     print("~~~")
#     print(processed, i)
#     for i in current_descendent:
#         if i in processed:
#             continue
#         s = x[i]
#         b1 = RNA.bp_distance(s1, s)
#         b2 = RNA.bp_distance(s2, s)
#         en = round(fc.eval_structure(s), 2)
#         candidates.append((s, i, b1, b2, en))
#         # print (s, i, b1, b2, en)

#     candidates.sort(key=lambda x: x[-1])
#     for s, i, b1, b2, en in candidates:

#         # sw = 32
#         # if RNA.bp_distance(s, s0) <= bp_dist0:
#         #     continue

#         # sw = RNA.bp_distance(s,s1)/2.0
#         sw = int(sw)
#         if sw<1: sw=1

#         rest1a = pathfinder.pathfinder(sequence, s, s1, search_width=sw).max_en
#         rest1b = pathfinder.pathfinder(sequence, s1, s, search_width=sw).max_en

#         # sw = RNA.bp_distance(s,s1)/2.0
#         sw = int(sw)
#         if sw<1: sw=1

#         rest2a = pathfinder.pathfinder(sequence, s, s2, search_width=sw).max_en
#         rest2b = pathfinder.pathfinder(sequence, s2, s, search_width=sw).max_en
        
#         rest1 = min(rest1a, rest1b)
#         rest2 = min(rest2a, rest2b)

#         # print ("fp call", b1, rest1)
#         # print (f's1 = "{s}"')
#         # print (f's2 = "{s1}"')


#         # print ("fp call", s, s2)

#         if max(rest1, rest2) < init_se:
#             print(s, i, b1, b2, b1+b2-bp_dist_all, en, max(rest1, rest2), "sw=", sw)

#             # rest1 = pathfinder.pathfinder(sequence, s, s1, search_width=sw, verbose=True).max_en
#             # rest2 = pathfinder.pathfinder(sequence, s, s2, search_width=sw, verbose=True).max_en

#         # sm = detect_local_minimum(fc, s)
#         # print (sm)

#     processed = processed | current_descendent
#     # print (last_descendent, current_descendent, current_descendent-last_descendent)

# # print()
# print(find_ancestor(s1_closest[0]))



x = RNAxplorer.repellant_sampling(fc)

print (x)