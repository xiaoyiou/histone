# This module is used to visualize the genes in network
from sklearn.metrics import jaccard_similarity_score
id2sax = dict()
sax2id = dict()
ind = 0
min_size = 10


for sax, v in all.saxGrps.items():
    if len(v) <= min_size:
        continue
    if sax not in sax2id:
        id2sax[ind] = sax
        sax2id[sax] = ind
        ind+=1

edges = [(0,0)] * (ind**2)
tail = 0
for i in xrange(ind):
    print i
    for j in xrange(i+1, ind):
        A, B = id2sax[i], id2sax[j]
        if jaccard_similarity_score(A,B) >= 0.8:
            edges[tail] = (i,j)
            tail += 1

