
import sax_tools as st
import cluster as cl
import pandas as pd
import saxnet as sn
import numpy as np
import gonet

M = 10
L = 21
alphabets = ['abcde'] * M
n_parts = [10] * M
levels = [3] * M
left = cl.createDataFrame(gData, L = L)
raw_data = pd.merge(left, young.labels, how='inner'\
                    , left_index=True, right_index=True)

all = cl.Spectrum(raw_data,mod_names)

all.dwt_transform(level=2)
#all.pp_transform(4,mode='paa')
all.perform_sax('abc')
all.ap_clustering()
all.getSaxGrps(thres=1)



"""
print "Clustering Finished for all genes"

lbs = young.labels
stress = all.sax_data.ix[lbs[lbs.stress == 1].index].dropna()
flower = all.sax_data.ix[lbs[lbs.flowerN == 1].index].dropna()
develop =  all.sax_data.ix[lbs[lbs.develop == 1].index].dropna()

sGrps = stress.groupby(by=stress.columns.tolist()).groups
fGrps = flower.groupby(by=flower.columns.tolist()).groups
dGrps = develop.groupby(by=flower.columns.tolist()).groups

temp1 = sorted(sGrps.items(), key = lambda x: len(x[1]))
temp2 = sorted(fGrps.items(), key = lambda x: len(x[1]))
temp3 = sorted(dGrps.items(), key = lambda x: len(x[1]))


for i in xrange(10):
    print "top %d" % i
    print "stress", temp1[-i][0]
    print "flower", temp2[-i][0]
    print "develop", temp3[-i][0]
    print "==================="

sRanks = dict()
for i, (k, v) in enumerate(temp1[::-1]):
    sRanks[k] = i + 1

fRanks = dict()
for i, (k, v) in enumerate(temp2[::-1]):
    fRanks[k] = i + 1
dRanks = dict()
for i, (k, v) in enumerate(temp3[::-1]):
    dRanks[k] = i + 1

patterns = list(set(sGrps.keys()) | set(fGrps.keys()) | set(dGrps.keys()))
rStress = []
rFlower = []
rDevelop = []
for p in patterns:
    rStress.append(sRanks[p] if p in sRanks else np.nan )
    rFlower.append(fRanks[p] if p in fRanks else np.nan)
    rDevelop.append(dRanks[p] if p in dRanks else np.nan)

temp = pd.DataFrame({'stress':rStress,'develop':rDevelop,'flower':rFlower}, index=patterns)
temp.to_csv('/Users/xiaoyiou/OneDrive/result/Feb_4/ranks.csv')
#snet.calcDists()
"""
#snet.exportNet('/Users/xiaoyiou/OneDrive/result/Feb_4/mod0.net')