# this module should take care of all weird stuff
from yxtools import getTabData
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sax_tools as st
import seaborn as sns
import sys

def getPatternInfo(data, path, obj, col_inds=[1, 2, 5, 6, 7], k=10, opath=None):
    """

    Parameters
    ----------
    data: should be the sorted saxGrps
    path: should be the path of gene names
    obj: the object of cdata.MData

    Returns: a list of patterns with sorted genes based on function labels
    -------

    """
    stdout = sys.stdout
    if opath:
        sys.stdout = open(opath, 'wb')
    names = dict()
    labels = obj.labels
    cols = labels.columns.tolist()
    col_names = [obj.classes[x] for x in col_inds]
    names = getTabData(path, delim='\t')
    res = []
    for pattern, glst in data:
        if len(glst) < k:
            break
        print '#'*50
        print pattern
        print '#' * 50 + '\n'
        t = labels.ix[glst]
        grps = t.groupby(by=col_names).groups
        for group, lst in grps.items():
            name = []
            for i,v in enumerate(group):
                if v != 0:
                    name.append(col_names[i])
            print '$' * 50
            print 'Labels: %s'%(' '.join(name) if len(name) else 'NO')
            for gene in lst:
                print gene, names.get(gene)
            print '$' * 50
            print ""
        print ''
    if opath:
        sys.stdout.close()
        sys.stdout = stdout


def getDistri(data, mod_names, length=21, nmod=10, window=1, label=None, norm=False, mode='average'):
    # transform the raw data into a list panda data frames where each column is a modification and each data frame
    # contains the data for a particular position (or window if window > 1)
    # if label is not None do filtering first
    temp = data.copy()
    res = []

    if label is not None :
        for l in label:
            temp = temp[temp[l] == 1]

    for i in xrange(length - window + 1):
        frame = {mod_names[x]: None for x in xrange(nmod)}
        for j in xrange(nmod):
            entry = temp[temp.modNum == j].ix[:, range(i, i+window)]
            if norm:
                entry = st.znormalization(entry)
            if mode == 'average':
                partial = entry.mean(axis=1)
            elif mode == 'max':
                partial = entry.max(axis=1)
            elif mode == 'all':
                # this one doesn't use combine the window it just concatenate the columns:
                partial = entry
            frame[mod_names[j]] = partial.values.flatten().tolist()
        res.append(pd.DataFrame(frame))

    return res


def visDists(A, method='pearson', size=None):
    # generate a bunch of plots for position-wise correlation
    for data in A:
        corrs = data.corr(method=method)
        plt.figure(figsize=size)
        temp = sns.heatmap(corrs)
        temp.set_xticklabels(temp.get_xticklabels(), rotation=30)
        temp.set_yticklabels(temp.get_yticklabels(), rotation=0)
        plt.show()











