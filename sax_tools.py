import numpy as np
import pandas as pd
import pywt
from scipy.stats import norm


# This module is mostly about transformation of data matrix
# it has no knowledge about genes names. It only know the matrix index

def znormalization(ts):
    """
    ts - each row is a histone modification record
    """
    mus = ts.mean(axis = 0)
    stds = ts.std(axis = 0)
    return (ts - mus) / stds

def paa_transform(ts, n_pieces):
    """
    ts: the columns of which are time series represented by e.g. np.array
    n_pieces: M equally sized piecies into which the original ts is splitted
    """
    splitted = np.array_split(ts, n_pieces, axis=1) ## along columns as we want
    return np.asarray(map(lambda xs: xs.mean(axis=1), splitted)).T

def pma_transform(ts, n_pieces):
    splitted = np.array_split(ts, n_pieces, axis=1)  ## along columns as we want
    return np.asarray(map(lambda xs: xs.max(axis=1), splitted)).T


def dist(dist, s1, s2 , p=2):
    assert len(s1) == len(s2)
    total = 0
    for i in xrange(len(s1)):
        total += dist.loc[s1[i], s2[i]] ** p
    return total ** (1.0/p)


def dwt(data, level=2, wavelet='db1'):
    res = []
    for _, row in enumerate(data):
        w = pywt.Wavelet(wavelet)
        coeff = pywt.wavedec(row, w, mode='constant')
        res.append(pywt.waverec(coeff[:level], w))
    return np.array(res,dtype=np.float)


def sax_transform(ts, alphabet, thrholds = None):
    """
    Before using this, the ts should be compressed (suing wavelet or paa)
    and then z transformed to ensure equal frequency assignments
    Return: (The sax representations matrix, distance matrix)
    """
    alphabet_sz = len(alphabet)
    if thrholds is None:
        thrholds = norm.ppf(np.linspace(1./alphabet_sz,
                                    1-1./alphabet_sz,
                                   alphabet_sz-1))

    def translate(ts_values):
        return np.asarray([(alphabet[0] if ts_value < thrholds[0]
                else (alphabet[-1] if ts_value > thrholds[-1]
                      else alphabet[np.where(thrholds <= ts_value)[0][-1]+1]))
                           for ts_value in ts_values])
    dist = np.zeros((alphabet_sz, alphabet_sz))
    for i in xrange(alphabet_sz - 1):
        for j in xrange(i + 1, alphabet_sz):
            if abs(i - j) <= 1:
                continue
            dist[i,j] = thrholds[j-1] - thrholds[i]
            dist[j,i] = dist[i,j]

    return (np.apply_along_axis(translate, 0, ts),\
            pd.DataFrame(dist, index=list(alphabet), columns=list(alphabet)),\
            thrholds)