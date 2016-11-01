# This script is to test the expanding idea
import vis
import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt
import mdata
import analysis as ana
from scipy import interp
from sklearn.datasets import make_classification
from sklearn.model_selection import StratifiedKFold as KFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve
from sklearn.naive_bayes import GaussianNB


X= binary[(binary.T!=0).any()]
#X=binary
X =X.ix[(X.index.str.contains('ATMG')==False)&(X.index.str.contains('ATCG')==False)]
glst_paths = ['defense.glst','development.glst'
              ,'flower.glst','flowering.glst','stress.glst'\
              ,'stimulus.glst','floweringN.glst']

inds = [0,1,0,0,1,1,0]

col_names = ['all','defense','develop','flowerB'\
             ,'flowerS','stress','stimulus','flowerN']
prefix ='genelst/'

dfLst = [X]+[ana.selectDataGL\
         (X,prefix+glst_paths[i],inds[i])\
         for i in range(len(glst_paths))]


gLsts = [x.index.tolist() for x in dfLst]

# -------------Actual Training and testing--------------
tprs = []
base_fpr = np.linspace(0, 1, 201)
fprs =[]
lw = 2
iters = 5
young =None
name = 'flowerN'

young = mdata.Mdata(X,mod_names)
for i in range(len(col_names)):
    young.addLabels(gLsts[i],col_names[i])
print "First iteration"
young.findPatterns(2)
print "calculating scores"
young.findPScores(mode='ar',ratios=False)
m = young.pgMap.copy()
XX = X
XX = m.apply(lambda row: row*young.ratios[name],axis=1).as_matrix()
yy = young.labels[name]
model = LogisticRegression(C=1e5).fit(XX, yy)
y_score = model.predict_proba(XX)
fpr, tpr, _ = roc_curve(yy, y_score[:, 1])
tpr = interp(base_fpr, fpr, tpr)
tpr[0] = 0.0

#for i in range(iters):
    







tprs.append(tpr)
    
    
tprs = np.array(tprs)



plt.plot(base_fpr, tpr, 'b--',label='%s, AUC: %.2f'%('some',metrics.auc(base_fpr,tpr)))
plt.plot([0, 1], [0, 1],'r--')
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 1.01])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.axes().set_aspect('equal', 'datalim')
plt.legend(loc="lower right")
plt.show()

    



