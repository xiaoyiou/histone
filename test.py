# This script is a test for hetero association learning using non-zero labels

from eonetwork import mapper,eval
import numpy as np
import matplotlib.pyplot as plt
import nengo
# from nengo import spa
#
#
# dim = 32
# vocab = spa.Vocabulary(dimensions=dim)
#
# words = ['ORANGE', 'APRICOT', 'CHERRY', 'STRAWBERRY', 'APPLE']
#
# for word in words:
#     vocab.parse(word)
#
#
#
#
#

res = {'pmloss':[], 'hloss':[], 'mlsetting':[], 'f1':[]}
nh_res = {'pmloss':[], 'hloss':[], 'mlsetting':[], 'f1':[]}

learner = mapper.Assoc(hopfield=True, keepz=False)
nh_learner = mapper.Assoc(hopfield=False, keepz=False)

A, B = learner.loaddata(all, go, 'F',useslim=False)

# removing all-zero entries
#inds = [x[0] for x in(B==1).any(axis=1).tolist()]
sizes = []

for i in xrange(10):
    inds = [x[0] for x in ((B == 1).sum(axis=1) > i).tolist()]
    if sum(inds) == 0:
        break
    a = A[inds, :]
    b = B[inds, :]
    sizes.append(a.shape[0])
    learner.train(a, b)
    nh_learner.train(a, b)
    zz = learner.recall(a)
    nh_zz = nh_learner.recall(a)
    for method in res:
        res[method].append(eval.eval(zz, b,type=method))
        nh_res[method].append(eval.eval(nh_zz, b,type=method))

