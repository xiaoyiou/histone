import cluster as cl
import pandas as pd
import datetime
from fastdtw import fastdtw
from functools import partial
from scipy.spatial.distance import euclidean

# Create the dataset first by combining labels with histone modifications


reload(cl)

L = 21
left = cl.createDataFrame(gData, L = L)
raw_data = pd.merge(left, young.labels, how='inner'\
                    , left_index=True, right_index=True)



all =  cl.Group(raw_data, 'all', mod_names, L, w=0)
# stress = cl.Group(raw_data, 'stress', mod_names, L, w=0)
#
#
#
#
# for i in xrange(stress.N):
#     test = stress.hmdata[i]
#     test.dwt()
# stress.cluster(2,5)





#stress.cluster(2,5)
#
#flower = cl.Group(raw_data, 'flowerN', mod_names, L, w=0)
#flowerData = gnet.GNet(gData,L, 10, flower.genes)
#flowerData.calcDists(euclidean)
#res = flowerData.getEdges()
# flower.cluster(2, 5)
#
# develop = cl.Group(raw_data, 'develop', mod_names, L, w=2)
# develop.cluster(2, 5)
#
# flowerPS = cl.findPatterns(flower)
# stressPS  = cl.findPatterns(stress)
# developPS = cl.findPatterns(develop)


# ------- pure randomness ------
"""
for i in xrange(5):
    temp = cl.Group(None, 'random_'+str(i),mod_names,L, w=1)
    temp.all_random(gData,300,mod_names)
    temp.cluster(2,10)
    temp.vis_all(mod_names)
    gc.collect()
"""

#Z = cl.HMData('develop', 3, L, data=raw_data,w=1)https://www.youtube.com/watch?v=dXbwC7Mdh6Qhttps://www.youtube.com/watch?v=dXbwC7Mdh6Qhttps://www.youtube.com/watch?v=dXbwC7Mdh6Q
#Z.opt_clustering((2,10))
# Only focusing on flowering for now
#ts_data = raw_data[((raw_data['stress'] == 1) | (raw_data['flowerN'] == 1))&(raw_data['modNum'] == 4)][range(L)].as_matrix()
#ts_data = raw_data[raw_data['modNum'] == 4][range(L)].as_matrix()

"""
C = tsc.ts_cluster(5)
C.k_means_clust(ts_data, 5, 5, progress=True)

clusters = C.assignments
Z = cl.getClusteringRes(clusters, ts_data)
plt.figure()
sns.tsplot(data=Z, time='pos', value='enrich', condition='cluster', unit='gene', ci=100)
"""