import matplotlib.pyplot as plt
import numpy as np
from scipy import interp

from sklearn.datasets import make_classification
from sklearn.model_selection import StratifiedKFold as KFold
#from sklearn.cross_validation import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve
from sklearn.naive_bayes import GaussianNB


name = 'salt'
method = 'HPSS'

m = young.pgMap.copy()

XPL = X[young.labels[name]==1]
# TT= XPL
TT =X.groupby(X.columns.tolist()).groups
TT =pd.DataFrame(sorted(TT,key=lambda k: len(TT[k]),reverse=True)[:10],columns=X.columns.tolist())




L = X.copy()
L['key']=L.index
mergedL = L.merge(TT,how='left',indicator=True)
mergedL = mergedL[mergedL['_merge']=='left_only']
R = XPL.copy()
R['key']=R.index
mergedR = R.merge(TT,how='left',indicator=True)
mergedR = mergedR[mergedR['_merge']=='left_only']




vindx =set(mergedL['key']).union(mergedR['key'])
#vindx=X.index


#XX =X.ix[vindx].as_matrix()
#XX = m.ix[vindx].apply(lambda row: row*young.ratios[name],axis=1).as_matrix()
XX = m.ix[vindx].apply(lambda row: row*young.pScore[name],axis=1).as_matrix()  
#temp =X[young.labels[name]==1]
#temp= temp.sum()/X.shape[0]
#XX=X.ix[vindx].apply(lambda row: row*temp,axis=1).as_matrix()




yy = young.labels[name].ix[vindx].as_matrix()



tprs = []
base_fpr = np.linspace(0, 1, 201)
fprs =[]
lw = 2

plt.figure(figsize=(10, 10))
idx = np.arange(0,len(yy))
for j in np.random.randint(0, high=10000, size=10):
    np.random.shuffle(idx)
    cv = KFold(n_splits=5,random_state=j)
    for i, (train, test) in enumerate(cv.split(XX,yy)):
#       model = LogisticRegression(C=1e5).fit(XX[idx][train], yy[idx][train])
        model = GaussianNB()
        model.fit(XX[idx][train],yy[idx][train])
        y_score = model.predict_proba(XX[idx][test])
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


plt.plot(base_fpr, mean_tprs, 'b--',label='%s, AUC: %.2f'%(method,metrics.auc(base_fpr,mean_tprs)))
plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='red', alpha=0.3)

plt.plot([0, 1], [0, 1],'r--')
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 1.01])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.axes().set_aspect('equal', 'datalim')
plt.legend(loc="lower right")
plt.show()
