# This script should be used to test the ratio changes
# between different gene lists
# we can "execfile" run this script after loadFiles.py

import analysis as ana
import rule as ru
from yxtools import dictInsert
from mdata import Mdata

glst_paths = ['defense.glst','development.glst'
              ,'flower.glst','flowering.glst','stress.glst'\
              ,'stimulus.glst']

inds = [0,1,0,0,1,1]

col_names = ['all','defense','develop','flowerB'\
             ,'flowerS','stress','stimulus']
prefix ='genelst/'

dfLst = [binary]+[ana.selectDataGL\
         (binary,prefix+glst_paths[i],inds[i])\
         for i in range(len(glst_paths))]


gLsts = [x.index.tolist() for x in dfLst]



train = True
if train:
        young =  mdata.Mdata(binary,mod_names)
        for i in range(len(col_names)):
                young.addLabels(gLsts[i],col_names[i])
        print "Searching for frequent itemsets"
        young.findPatterns(2)
        print "calculating scores"
        young.findPScores()

print "predicting"
result = young.predictALL(mode='norm')
plt.clf()

for name in young.classes:
    scores =result[name]
    y = young.labels[name]
    fpr,tpr,thresholds=metrics.roc_curve(y,scores)
    plt.plot(fpr, tpr, label='ROC curve for %s'%(name))
    

plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="upper left")
plt.show()

