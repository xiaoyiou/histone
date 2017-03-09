# This script is to test the performance of classification using new feature space
# this should be executed after "groupsax_testing,py" and "annotation_test"
from eonetwork import mapper
import mdata
import pandas as pd
method = 'SAXGRP'
name = 'flowerN'
temp = mapper.Assoc()
A, B = temp.loaddata(all, go, 'F')

A[A==-1] = 0




"""
plt.figure(figsize=(10, 10))
idx = np.arange(0, len(yy))
for j in np.random.randint(0, high=10000, size=10):
    np.random.shuffle(idx)
    cv = KFold(n_splits=5,random_state=j)
    for i, (train, test) in enumerate(cv.split(A,yy)):
        model = LogisticRegression(C=1e2).fit(A[idx], yy[idx])
#        model = GaussianNB()
        model.fit(A[idx], yy[idx])
        y_score = model.predict_proba(A[idx][test])
        fpr, tpr, _ = roc_curve(yy[idx][test], y_score[:, 1])

#        plt.plot(fpr, tpr, 'b', alpha=0.05)
        tpr = interp(base_fpr, fpr, tpr)
        tpr[0] = 0.0
        tprs.append(tpr)

tprs = np.array(tprs)
mean_tprs = tprs.mean(axis=0)
std = tprs.std(axis=0)

tprs_upper = np.minimum(mean_tprs + std, 1)
tprs_lower = mean_tprs - std
plt.plot(base_fpr, mean_tprs, 'b--', label='%s, AUC: %.2f' % (method, auc(base_fpr, mean_tprs)))
plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='red', alpha=0.3)

plt.plot([0, 1], [0, 1], 'r--')
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 1.01])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.axes().set_aspect('equal', 'datalim')
plt.legend(loc="lower right")
plt.show()
"""