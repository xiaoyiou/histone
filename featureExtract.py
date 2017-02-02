"""
Created on Tue Nov 29 13:04:46 2016

@author: xiaoyiou
"""
import pandas as pd
import numpy as np
import copy
from fastdtw import fastdtw
from numpy.linalg import norm
from scipy.spatial.distance import euclidean

# This module deals with feature extraction from raw enrichment data


def createData(gData):
    res = {}
    for gene, mod in gData:
        if not gene in res:
            res[gene] = {}
        res[gene][mod] = gData[(gene,mod)]
    return res


def binarize(data, threshs):

    res = copy.deepcopy(data)
    for key in data:
        for mod in data[key]:
            res[key][mod] = [ 1 if x > threshs[mod] else 0 for x in res[key][mod]]
    return res


def accumuPeaks(bData):
    res = copy.deepcopy(bData)
    for key in bData:
        for mod in bData[key]:
            res[key][mod] = [sum(bData[key][mod][:i+1]) for i in xrange(len(bData[key][mod]))]
    return res


def createDistMatrix(data):
    N = len(data)
    res = np.zeros([N,N])
    names = {}
    memo = {}
    i = 0
    for gene in data:
        names[i] = gene
        i +=1
    for i in xrange(N):
        for j in xrange(0,N):
            if i == j:
                res[i][j] = 0
            else:
                total = 0
                L = len(data[names[i]])
                dist = 0
                for k in xrange(L):
                    x = data[names[i]][k]
                    y = data[names[j]][k]
                    xkey = ''.join(str(c) for c in x)
                    ykey = ''.join(str(c) for c in y)
                    key1 = xkey+'+'+ykey
                    key2 = ykey+'+'+xkey
                    if key1 in memo:
                        dist = memo[key1]
                    elif key2 in memo:
                        dist = memo[key2]
                    else:
                        dist,_ = fastdtw(x,y, dist= euclidean)
                        memo[key1] = dist
                        memo[key2] = dist
                    total += dist
                res[i][j] = total/float(L)
    return (res,names)

