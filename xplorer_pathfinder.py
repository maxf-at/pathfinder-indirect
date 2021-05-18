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


# s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..' #s1
# s2       = '....(((.((((...((((....))))(((((.(((......)))))))))))))))......((((((((((........).)))))))))........' #mfe

# s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..' #s1
# s2       = '..........(((...((((((..(((((((..(((......))).)))..((((...((....(((.....)))....)))))))))))))))))))..' #A

# s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..' #s1
# s2       = '..........(((...((((((..((((..((.(((......))).))....(((...((....(((.....)))....))))).)))))))))))))..' #B

# s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..' #s1
# s2       = '..........(((...((((((..((((.(((.(((......))).)))(((..((..............))..)))........)))))))))))))..' #C
# s1       = '..........(((...((((((..((((.(((.(((......))).)))(((..((..............))..)))........)))))))))))))..' #C
# s2       = '..........(((...((((((..((((..((.(((......))).))....(((...((....(((.....)))....))))).)))))))))))))..' #B
# s1       = '..........(((...((((((..((((..((.(((......))).))....(((...((....(((.....)))....))))).)))))))))))))..' #B
# s2       = '....(((.((((...((((....))))(((((.(((......)))))))))))))))......((((((((((........).)))))))))........' #mfe



# x = callRNAxplorer(sequence, s1, s2, 2)

# /**
#  * Compute a distance based path between structure s1 and s2.
#  * @param seq - the RNA sequence (ACGU)
#  * @param s1 - first structure in dot-bracket format.
#  * @param s2 - second structure in dot-bracket format.
#  * @param iteration - the current iteration in a series of recursive calls. It should be initialized with 0.
#  * @param maxIterations - maximal iterations of structurewalks.
#  * @param maxKeep - how many structures are being kept (get_path in findpath.c)
#  * @param maxStorage - maximal number of meshpoints (see insert_meshpoint in meshpoint.c)
#  */


x = RNAxplorer.levelSaddlePoint(sequence, s1, s2, 100, 5000, 1, 100)


# x = RNAxplorer.levelSaddlePoint2(sequence, s1, s2, 2, 500, 100, 100, 20, 10)

print (x)