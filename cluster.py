"""
This module is used to do clustering of genes based on modification levels

using the gene list and modification pattern to filter out the genes

"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsClassifier as KNN
from sklearn.metrics import silhouette_score
from fastdtw import fastdtw
from scipy.spatial.distance import euclidean
from sklearn.cluster import KMeans, DBSCAN
from scipy import signal
import peakutils
import gc
from functools import partial
from multiprocessing import Process, Manager
import saxnet as sn
import sax_tools as st


class HMData(object):
    """
    HMData is an encapsulation of the hitone modification enrichment levels
    of a particular mod and label
    """
    def __init__(self, label, modNum, L, data=None, w=3):
        self.label = label
        self.mod = modNum
        self.L = L
        self.data = None
        self.genes = None
        self.N = 0
        self.clusters = None
        self.centroids = None
        self.k = 0
        self.w = w
        self.sims = None
        self.score = -1
        self.ass = None
        self.distr = None
        self.learner = None
        self.dists = None
        self.plateau = None
        self.peaklocs = dict()
        self.approx = None
        if data is not None:
            self.extractData(data)


    def getRandomData(self, gData, L, size, mod):
        self.genes = []
        data = []
        allgs = []
        for gene, _ in gData.keys():
            allgs.append(gene)
        allgs = list(set(allgs))

        candidates = np.random.choice(allgs, replace=False, size=size)
        for gene in candidates:
            if len(gData[(gene,mod)]) != L:
                continue
            self.genes.append(gene)
            data.append(gData[(gene,mod)])
        self.data =pd.DataFrame(data).as_matrix()
        self.N = self.data.shape[0]
        self.__smooth()

    def calSims(self):
        # calculate the similarity matrix
        self.sims = np.zeros((self.N, self.N))
        C = tsc.ts_cluster(3)
        for i in xrange(self.N):
            for j in xrange(i+1, self.N):
                self.sims[i, j]=C.DTWDistance(self.data[i, :],\
                                             self.data[j, :], w=self.w)
                self.sims[j, i] = self.sims[i, j]

    def __smooth(self, window):
        if self.w == 0:
            return

        def __convol(row):
            w = self.w
            window = np.ones(w) / float(w)
            return np.convolve(row, window, 'same')
        self.data = np.apply_along_axis(__convol, 0, self.data)


    def extractData(self, raw_data):
        temp = raw_data[(raw_data[self.label] == 1) & (raw_data['modNum'] == self.mod)]
        self.data = temp[range(self.L)].as_matrix()
        self.genes = temp.index.tolist()
        self.N = len(self.genes)
        #self.__smooth()
        print "loading data finished"

    def get2DMax(self):
        """
        for each gene, calculates the loc of the peak and the amplitude of the peak
        -------
        """
        pass

    def __cluster(self, k):

        #C = tsc.ts_cluster(k)
        #C.k_means_clust(self.data, n_iter, window, progress=True)

        kcl = KMeans(n_clusters=k)
        print "starting to fit k:" + str(k)
        if self.approx is not None:
            kcl.fit(self.approx)
        else:
            kcl.fit(self.data)
        print "clustering (k=%d) finished" % (k)
        return kcl


    def dbcluster(self, **args):
        dbcl = DBSCAN(**args)
        if self.approx is not None:
            dbcl.fit(self.approx)
        else:
            dbcl.fit(self.data)
        return dbcl

    def __getScore(self, C):
        """

        Parameters
        ----------
        C: the ts_cluster object

        Returns
        -------

        """
        res = 0
        labels = C.labels_
        if self.approx is not None:
            res = silhouette_score(self.approx, labels, metric='euclidean')
        else:
            res = silhouette_score(self.data, labels, metric='euclidean')
        return res

    def opt_clustering(self, krange):
        low, high = krange

        for k in xrange(low,high + 1):
            C = self.__cluster(k)
            score = self.__getScore(C)
            if score > self.score:
                self.k = k
                self.score = score
                self.centroids = C.cluster_centers_
                self.clusters = C.labels_
                self.learner = C
            gc.collect()
        print "best k=%d for [%s,%d]" % (self.k, self.label, self.mod)
        self.getAssignment()
        self.getDistr()

    def getAssignment(self):
        """

        Returns Nothing, only calculate the assignments for all the genes using
        {1: [AT1G0001, AT2G2002. etc] }
        -------

        """
        self.ass = {}
        for ind, c in enumerate(self.clusters):
            if c not in self.ass:
                self.ass[c] = []
            self.ass[c].append(self.genes[ind])

    def getDistr(self):
        self.distr = [0] * len(self.ass)
        for i in self.ass:
            self.distr[i] = len(self.ass[i]) / float(self.N)

    def vis_cluster(self, mod_names):

        cluster = {}
        for ind, v in enumerate(self.clusters):
            if v not in cluster:
                cluster[v] = []
            cluster[v].append(ind)
        pData = None
        if self.approx is None:
            pData = getClusteringRes(cluster, self.data)
        else:
            pData = getClusteringRes(cluster, self.approx)
        plt.figure()
        sns.tsplot(data=pData, time='pos', value='enrich',\
                   condition='cluster', unit='gene', ci=95.0)
        plt.title("Function: %s, mod: %s" % (self.label,mod_names[self.mod]))
        plt.show()

    def vis_centroids(self, mod_names):

        plt.figure()
        for ind, v in enumerate(self.centroids):
            plt.plot(self.centroids[ind], label='Cluster %d' % (ind), linewidth=10*self.distr[ind])
        plt.legend()
        plt.show()

    def peakify(self, thres=0.5):
        self.plateau = np.zeros((self.N, self.L))
        for ind, v in enumerate(self.data):
            inds = set(peakutils.indexes(v,thres=thres, min_dist=self.w))
            for j in xrange(self.L):
                if j not in inds:
                    self.data[ind,j] = 0
                    self.peaklocs[ind] = inds
            for i in inds:
                for j in xrange(max(0, i - self.w), min(self.L, i+self.w+1)):
                    self.plateau[ind, j] = self.data[ind, i]

    def calcAllDists(self, func, n_jobs =4):
        self.dists = dict()
        def foo(res, data, low, high, func):
            for x in xrange(low,high):
                v = data[x, :]
                for y,u in enumerate(self.data):
                    edge = frozenset({x, y})
                    if edge in res or x == y:
                        continue
                    res[edge]= func(v,u)

        temp = Manager().dict()
        step = self.N / n_jobs
        lows = [x*step for x in xrange(n_jobs)]
        highs = [min(x+step,self.N) for x in lows]
        workers = []
        for i in range(n_jobs):
            p = Process(target=foo, args=(temp, self.data,lows[i],highs[i],func))
            p.start()
            p.join()
            workers.append(p)
        for key in temp.keys():
            self.dists[key] = temp[key]
        del temp
        for p in workers:
            if p.is_alive():
                p.terminate()

    def calcAllDistsDTW(self, power=1, n_jobs =4):
        self.dists = dict()
        def foo(res, low, high):
            for x in xrange(low,high):
                inds = self.peaklocs[x]
                v = self.data[x,:]
                for y,u in enumerate(self.plateau):
                    edge = frozenset({x, y})
                    if edge in res or x == y:
                        continue
                    res[edge] = sum([min(v[i], (u[i]-v[i]))**power for \
                                      i in inds])**(1.0/power)

        temp = Manager().dict()
        step = self.N / n_jobs
        lows = [x*step for x in xrange(n_jobs)]
        highs = [min(x+step,self.N) for x in lows]
        for i in range(n_jobs):
            p = Process(target=foo, args=(temp,lows[i],highs[i]))
            p.start()
            p.join()
        for key in temp.keys():
            self.dists[key] = temp[key]
        del temp
        # for p in workers:
        #     if p.is_alive():
        #         p.terminate()

class Group(object):

    def __init__(self, raw_data, func_name, mod_names, L, w=3):
        self.name = func_name
        self.hmdata = []  # A list of HMData objects indexed by mod_numbers
        self.mod_names = list(mod_names)
        self.w = w
        self.L = L
        self.genes = None
        self.distrs = []
        self.N = len(mod_names)
        if raw_data is None:
            return
        for i in xrange(len(mod_names)):
            temp = HMData(func_name,i,L, data=raw_data, w=self.w)
            if self.genes is None:
                self.genes = list(temp.genes)
            self.hmdata.append(temp)
            print "loaded data for %s" % (mod_names[i])

    def cluster(self, klow, khigh):
        for obj in self.hmdata:
            obj.opt_clustering((klow, khigh))
            print "clustered data for %s" %(self.mod_names[obj.mod])
            self.distrs.append(obj.distr)

    def all_random(self, gData, size, mod_names):

        for i in xrange(len(mod_names)):
            temp = HMData('random',i, self.L,w=self.w)
            temp.getRandomData(gData, self.L, size, i)
            self.hmdata.append(temp)
        gc.collect()


    def vis_all(self,mod_names, nc = 5):
        f, axarr = plt.subplots(self.N/nc + 0 if self.N % nc ==0 else 1, nc)
        for i in xrange(self.N):
            for ind, v in enumerate(self.hmdata[i].centroids):
                axarr[i/nc, i % nc].plot(v, label='Cluster %d' % ind, \
                linewidth=10*self.hmdata[i].distr[ind])
                axarr[i/nc, i % nc].set_title(mod_names[i] \
                                              + '_%.2f'%self.hmdata[i].score)
                axarr[i/nc, i % nc].legend()
        plt.show()


    def getCatData(self):
        """
        After the clustering, using the assigned labels to discretize the data into
        categorical dataset
        Returns
        -------
        """
        data = pd.DataFrame(np.zeros((len(self.genes),self.N)), index = self.genes)

        for i in xrange(self.N):
            ass = self.hmdata[i].ass
            for key in ass:
                for ind in ass[key]:
                    data.loc[ind, i] = key
        return data

class Spectrum(object):
    def __init__(self, raw_data, mod_names, L=21):
        """

        Parameters
        ----------
        rawData:        the rawdata
        mod_names:      modification English names
        L:              the number of data points in each time series
        alphabets:      number of vertical lines
        """
        self.saxnets = []  # A list of HMData objects indexed by mod_numbers
        self.mod_names = list(mod_names)
        self.L = L
        self.genes = None
        self.M = len(mod_names)
        self.N = 0
        self.saxGrps = None
        self.saxNet = None

        if raw_data is None:
            return
        for i in xrange(len(mod_names)):
            temp = HMData('all', i, L, data=raw_data, w=0)
            if self.genes is None:
                self.genes = list(temp.genes)
                self.N = len(self.genes)
            self.saxnets.append(sn.SaxNet(temp.data))
            print "loaded data for %s" % (mod_names[i])

    def dwt_transform(self, level=2):
        for i in xrange(self.M):
            self.saxnets[i].dwt_transform(level=level)

    def pp_transform(self, n_parts, mode='paa'):
        for i in xrange(self.M):
            self.saxnets[i].pp_transform(n_parts, mode=mode)

    def perform_sax(self, alphabets):
        for i in xrange(self.M):
            self.saxnets[i].perform_sax(alphabets)

    def ap_clustering(self):
        for i in xrange(self.M):
            self.saxnets[i].ap_clustering()
            print "clustering done for mod %d" % i
    def calc_dists(self):
        for i in xrange(self.M):
            self.saxnets[i].calcDists()


    def getSaxGrps(self, thres = 5):
        all_saxes = np.zeros((self.N, self.M))
        for i in xrange(self.M):
            all_saxes[:, i] = self.saxnets[i].ass

        self.sax_data = pd.DataFrame(all_saxes, index=self.genes, columns=self.mod_names)
        self.saxGrps = self.sax_data.groupby(by=self.sax_data.columns.tolist()).groups
        for k,v in self.saxGrps.items():
            if len(v) < thres:
                del self.saxGrps[k]


    def vis_all(self,mod_names, nc = 5, size=None):
        f, axarr = plt.subplots(self.M/nc + 0 if self.N % nc ==0 else 1, nc)
        for i in xrange(self.M):
            distr = map(len, self.saxnets[i].clusters.values())
            print distr
            for ind, v in enumerate(self.saxnets[i].cluster_centers):
                axarr[i/nc, i % nc].plot(v, label='Cluster %d' % ind, \
                linewidth=10*distr[ind]/1.0/self.N)
                axarr[i/nc, i % nc].set_title(mod_names[i])
                axarr[i/nc, i % nc].legend()
        if size:
           f.set_size_inches(size[0],size[1])
        plt.show()




def flattenData(enrich,labels,pattern=None,method='original',w=0):
    """
    Parameters
    ----------
    enrich : the nested dict {gene:{mod:[...]}} of enrichment level
    pattern : A subset of modifications
    method : 'original' means we don't change the input data,
    'min' will take the minimum from the window ( windows size = w)
    'max' will take the maximum from the window
    'average" will take the average from the window
    Returns: The X and y
    -------
    """
    name_map = {}
    data = []
    ind = 0
    for gene in enrich:
        name_map[ind] = gene
        row = []
        if pattern:
            for p in pattern:
                row += __warping(enrich[gene][p], method=method, w=w)
        else:
            for p in enrich[gene]:
                row += __warping(enrich[gene][p], method=method, w=w)

        data.append(row)
        ind += 1

    return np.array(data), labels[labels.columns[1:]].as_matrix(), name_map


def createDataFrame(enrich, L = 21):
    """

    Parameters
    ----------
    enrich: the gData dict {(gene, modNum): []} dictionary of raw enrichment levels
    labels: the 1/0 labels DataFrame from cData object

    Returns a Pandas DataFrame with the following structure

    pos1, pos2, ... pos L , stress, flowering, ... development, mod_Num
    -------
    """
    data = {}
    for i in xrange(L):
        data[i] = []
    indexes = []
    data['modNum'] = []
    for gid, modNum in enrich:
        if len(enrich[(gid, modNum)]) != L:
            continue
        indexes.append(gid)
        data['modNum'].append(modNum)
        for i in xrange(L):
            data[i].append(enrich[(gid, modNum)][i])
    return pd.DataFrame(data, index = indexes)



def singleFlatten(enrich,labels, label, M= 10, L = 21,pattern=None, method = 'original', w = 0):
    name_map = {}
    data = []
    ind = 0
    indx = []
    width = 0
    y = []
    if pattern is None:
        width = M*L if method == 'original' else M * (L-w + 1)
    else:
        width = len(pattern) * L if method == 'original' else len(pattern)*(L-w + 1)

    for gene in enrich:
        name_map[ind] = gene
        row = []
        if pattern:
            for p in pattern:
                row += __warping(enrich[gene][p], method=method, w=w)
        else:
            for p in enrich[gene]:
                row += __warping(enrich[gene][p], method=method, w=w)
        if len(row) == width:
            data.append(row)
            indx.append(ind)
            y.append(labels.ix[gene,label])
            ind += 1
    return np.array(data), np.array(y), name_map



def getTopPattern(obj, label, k = 5):
    """

    Parameters
    ----------
    obj : the Mdata object
    label: the label we are interested in
    k: the top k patterns to consider

    Returns
    -------

    """
    ptns = obj.pScore.sort_values(by=label, ascending=False)[:k].index.tolist()
    res = [list(x) for x in ptns]
    res = sum(res,[])
    return list(set(res))



def trainKNN(X, y, k, weights='uniform', leaf_size =30, p=2):
    """

    Parameters
    ----------
    X: A list of lists to represent the enrichment levels
    y: A list of labels(the index of mods), -1 menas unlabelled
    k: n_neighbors
    weights: default is uniform. 'distance' means inverse of distance, or a callable function
              which accepts any array of distances and return an array of the same shape
    leaf_size: determines the size of brute-force nodes
    p: the dimension of distance p =1 means manhattan, p = 2 means euclidean distance

    Returns
    -------

    """
    if len(X) != len(y):
        print "Unmatching dimensions of input and labels"
        return

    clf = KNN(n_neighbors=k, weights=weights,leaf_size=leaf_size,p=p,n_jobs=8)
    clf.fit(X,y)
    return clf



def getClusteringRes(clusters, ts_data):
    n, m = ts_data.shape
    data = {'pos':[], 'gene':[],'cluster':[],'enrich':[]}
    for ind, v in enumerate(ts_data):
        for i in xrange(len(v)):
            data['pos'].append(i)
            data['enrich'].append(v[i])
            data['gene'].append(ind)
            data['cluster'].append(-1)
    res = pd.DataFrame(data)
    for key in clusters:
        for i in clusters[key]:
            res.loc[res.gene == i, ['cluster']] = key
    return res
    

def findPatterns(group):
    """

    Parameters
    ----------
    data

    Returns
    -------

    """

    data = group.getCatData()

    return sorted(data.groupby(data.columns.tolist()).groups.items(),\
                            key = lambda x: len(x[1]))

