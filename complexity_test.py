# This module is to test the complexity of data modeling
import cluster as cl
import pandas as pd
alphabets = 'abcde'
levels = [1,2,3,4]
alphas = [2,3,4,5]

ngrps = [[0]*4 for _ in xrange(len(levels))]


for i in xrange(len(levels)):
    for j in xrange(len(alphas)):
        test = cl.Spectrum(raw_data, mod_names)
        test.dwt_transform(level=levels[i])
        test.perform_sax(alphabets[:alphas[j]])
        test.ap_clustering()
        test.getSaxGrps(thres=5)
        ngrps[i][j] = len(test.saxGrps)