# This script should be used to test the ratio changes
# between different gene lists
# we can "execfile" run this script after loadFiles.py

import analysis as ana
import rule as ru
from yxtools import dictInsert
import mdata 
import matplotlib.pyplot as plt
from sklearn import metrics
import vis

X= binary[(binary.T!=0).any()]
#X=binary
X =X.ix[(X.index.str.contains('ATMG')==False)&(X.index.str.contains('ATCG')==False)]
glst_paths = ['defense.glst','development.glst'
              ,'flower.glst','flowering.glst','stress.glst'\
              ,'stimulus.glst','floweringN.glst','salt.glst','heat.glst']

inds = [0,1,0,0,1,1,0,0,0]

col_names = ['all','defense','develop','flowerB'\
             ,'flowerS','stress','stimulus','flowerN',   'salt','heat']
prefix ='genelst/'

dfLst = [X]+[ana.selectDataGL\
         (X,prefix+glst_paths[i],inds[i])\
         for i in range(len(glst_paths))]


gLsts = [x.index.tolist() for x in dfLst]


young = None
train = True
if train:

        young =  mdata.Mdata(X,mod_names)
        for i in range(len(col_names)):
                young.addLabels(gLsts[i],col_names[i])
        print "Searching for frequent itemsets"
        young.findPatterns(2)
        print "calculating scores"
        young.findPScores(mode='ar',ratios=False)

print "predicting"
#result = young.predictALL(young.binary,mode='norm')


"""
plt.clf()
for name in young.classes[1:]:
    scores =result[name]
    y = young.labels[name]
    fpr,tpr,thresholds=metrics.roc_curve(y,scores)
    plt.plot(fpr, tpr, label='%s, AUC: %.2f'%(name,metrics.auc(fpr,tpr)))
    

plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="upper left")
plt.show()

"""
