# This module should be used as a playground for enrichment clustering project

import featureExtract as fe
import cluster
import numpy as np
from sklearn.metrics import roc_curve
from sklearn.model_selection import StratifiedKFold as KFold
from scipy import interp
from sklearn import metrics

filter_pattern = None
#filter_pattern = cluster.getTopPattern(young,'flowerN',5)

def getRocs(k, X, y, n_splits=3, repeats=3):
    tprs = []
    base_fpr = np.linspace(0,1,201)
    fprs = []
    print X.shape, y.shape
    for _ in xrange(repeats):
        cv = KFold(n_splits=n_splits)
        for i, (train, test) in enumerate(cv.split(X, y)):
            model = cluster.trainKNN(X[train],y[train],k,weights='uniform',leaf_size=20,p=1)
            model.fit(X[train],y[train])
            probs = model.predict_proba(X[test])
            fpr, tpr, _ = roc_curve(y[test], probs[:,1])
            tpr = interp(base_fpr, fpr, tpr)
            tpr[0] = 0.0
            tprs.append(tpr)
    tprs = np.array(tprs)
    mean_tprs = tprs.mean(axis = 0)
    std = tprs.std(axis = 0)
    tprs_upper = np.minimum(mean_tprs + std, 1)
    tprs_lower = mean_tprs - std
    return mean_tprs, std, tprs_upper, tprs_lower, base_fpr



# loadFiles.py, cData should be run at this point
sub_gData = {}
for gene,mod in gData:
    if gene in young.labels.index:
        sub_gData[(gene,mod)] = gData[(gene,mod)]

enrich = fe.createData(sub_gData)
# only focus on flowering at the moment
cluster_x, cluster_y, name_map = cluster.singleFlatten(enrich, young.labels,\
                    'all', method='original', w=5, pattern=filter_pattern)
"""
# Now testing the original data set
mean_tprs, std_tprs, tprs_upper, tprs_lower, base_fpr =\
getRocs(100, cluster_x, cluster_y)

plt.figure()
plt.plot(base_fpr, mean_tprs, 'b',label='AUC: %.2f'%(metrics.auc(base_fpr,mean_tprs)))
plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='red', alpha=0.3)

plt.plot([0, 1], [0, 1],'r--')
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 1.01])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.axes().set_aspect('equal', 'datalim')
plt.legend(loc="lower right")
plt.show()
"""