"""
This module aims at creating the gene networks using histone modifications
"""

import numpy as np
from igraph import *
import gc

class GNet(object):

    def __init__(self, gData, L, M, filt=None):
        self.data = dict()
        self.L = L
        self.M = M
        for gene, mod in gData:
            if filter and gene not in filt:
                continue
            if len(gData[(gene,mod)]) != L:
                continue

            if gene not in self.data:
                self.data[gene] = np.zeros((M,L))
            self.data[gene][mod,:] = gData[(gene,mod)]
        self.dists = dict()

    def calcDists(self, func):
         """

         Parameters
         ----------
         func  Should be a partial function defined outsied the class which only takes two lists of the
        same length and calculate the distance

         Returns
         -------

         """
         for source in self.data:
             for target in self.data:
                edge = frozenset({source, target})
                if edge in self.dists or source == target:
                    continue
                self.dists[edge] = [0] * self.M
                for i in xrange(self.M):
                    self.dists[edge][i] = func(self.data[source][i,:]\
                                               ,self.data[target][i,:])
    def getEdges(self,alpha=1,selector=None):
        """
        If selector is None: then we use the average of distances
        if it's a number 0,1,2...M-1 then just use the corresponding distance
        Parameters
        ----------
        selector

        Returns
        -------
        """

        temp = self.dists.copy()
        if selector is None:
            for key in temp:
                temp[key] = np.mean(temp[key])
        else:
            for key in temp:
                temp[key] = temp[key][selector]
        thresh =np.mean(temp.values()) - alpha*np.std(temp.values())

        return [tuple(key) for key in temp if temp[key] <= thresh]





