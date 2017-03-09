import sax_tools as st
import pandas as pd
from sklearn.cluster import AffinityPropagation, KMeans, DBSCAN
import numpy as np

from scipy.spatial.distance import euclidean
class SaxNet(object):

    def __init__(self, data):
        self.data = data
        self.alphabets = None
        self.n_parts = None
        self.N, self.M = data.shape
        self.lookup = None   # the distance matrix for SAX method
        self.indexes = None  #A the SAX groups
        self.dists = None
        self.cdata = None # the compressed data
        self.bins = None
        self.avdata = dict()
        self.saxsize = dict()
        # The info about the clustering result
        self.n_clusters = 0
        self.clusters = None
        self.cluster_centers = None
        self.cluster_sax = None
        self.asax_data = None
        self.ass = None     # the cluster assignments for each data point, a long list of
                            # sax ids

    def pp_transform(self, n_parts, mode='paa', ztrans=True):
        self.n_parts = n_parts
        if mode == 'paa':
            self.cdata = st.paa_transform(self.data, self.n_parts)
        elif mode == 'pma':
            self.cdata = st.pma_transform(self.data, self.n_parts)
        if ztrans:
            self.cdata = st.znormalization(self.cdata)

    def dwt_transform(self, level=2, wave='haar', ztrans=True):
        self.cdata = st.dwt(self.data, level=level, wavelet=wave)
        if ztrans:
            self.cdata = st.znormalization(self.cdata)
        self.n_parts = self.cdata.shape[1]

    def perform_sax(self, alphabets):
        self.alphabets = alphabets
        if self.cdata is not None:
            sax_data, self.lookup, self.bins = st.sax_transform(self.cdata, self.alphabets)
        else:
            sax_data, self.lookup, self.bins = st.sax_transform(self.data, self.alphabets)
        C = self.cdata.shape[1] if self.cdata is not None else self.data.shape[1]
        self.indexes = pd.DataFrame(sax_data).groupby(by=range(C)).groups
        for sax in self.indexes:
            temp = self.cdata if self.cdata is not None else self.data
            self.saxsize[sax] = len(self.indexes[sax])
            self.avdata[sax] = temp[self.indexes[sax], :].mean(axis=0)

    def ap_clustering(self, preference=None, precomputed = True):
        if precomputed:
            cls = AffinityPropagation(preference=preference, affinity='precomputed')
            data = self.pw_dists()
        else:
            AffinityPropagation(preference=preference)
            data = self.avdata.values() if len(self.avdata) else \
                (self.cdata if self.cdata is not None else self.data)
        cls.fit(data)
        reps = self.indexes.keys()
        self.cluster_sax = [reps[i] for i in cls.cluster_centers_indices_]
        self.cluster_centers = [self.avdata[sax] for sax in self.cluster_sax]
        self.clusters = dict()
        for ind, label in enumerate(cls.labels_):
            sax = self.cluster_sax[label]
            if sax not in self.clusters:
                self.clusters[sax] = []
            self.clusters[sax] += self.indexes.values()[ind]
        self.asax_data = dict()
        for sax in self.clusters:
            temp = self.cdata if self.cdata is not None else self.data
            self.asax_data[sax] = temp[self.clusters[sax], :].mean(axis=0)
        self.ass = [0] * self.N
        for sax in self.cluster_sax:
            v = self.cluster_sax.index(sax)
            for ind in self.clusters[sax]:
                self.ass[ind] = v

    def pw_dists(self):
        N = len(self.indexes)
        dist = np.zeros((N,N))
        for i in xrange(N):
            for j in xrange(i+1, N):
                x, y = self.indexes.keys()[i], self.indexes.keys()[j]
                dist[i,j] = st.dist(self.lookup,x,y)
                dist[j,i] = dist[i,j]
        return dist




    def calcDists(self):
        self.dists = dict()
        for x in self.cluster_sax:
            for y in self.cluster_sax:
                if x == y:
                    self.dists[frozenset({x,y})] = 0
                key = frozenset({x,y})
                if key in self.dists:
                    continue

                self.dists[key] = st.dist(self.lookup, x, y)

    def getDist(self,n1, n2):
        x, y = self.cluster_sax[n1],self.cluster_sax[n2]
        return self.dists[frozenset({x,y})]



    def exportNet(self, path = None):
        if path is None:
            for a, b in self.dists:
                print "%s\t%s\t%.3f" % (''.join(a), ''.join(b), self.dists[frozenset({a, b})])
        else:
            with open(path,'wb') as f:
                f.write("*Vertices %d\n" % len(self.ninds))
                for i, v in enumerate(self.ninds):
                    f.write("%d\t%s\t%d\n" % (i + 1, v, self.nodes[i].N))
                f.write("*Arcs %d\n" % len(self.dists))
                for a, b in self.dists:
                    sim = 1/(1+self.dists[frozenset({a, b})])
                    if sim < 0.8:
                        continue
                    f.write( "%d\t%d\t%.3f\n" % (a + 1, b + 1, sim))
                f.close()
