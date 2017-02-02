import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
# at this point, the enrich_test should be run already

compressed_x = TSNE().fit_transform(cluster_x)
colors = [0] * compressed_x.shape[0]

cmap = {}
cmap[0] = 'w'
cmap[1] = 'c'
cmap[2] = 'g'
cmap[3] = 'r'
cmap[4] = 'r'
cmap[5] = 'k'
cmap[6] = 'y'
cmap[7] = 'r'

for i in xrange(compressed_x.shape[0]):
    row = young.labels.ix[name_map[i]]
    print row
    if sum(row) == 2:
        colors[i] = cmap[row.tolist()[1:].index(1)+1]
    else:
        colors[i] = cmap[0]

