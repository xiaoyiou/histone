# This script is to test
import vis
import rule as ru
import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt


if 'asocL' in locals():
    print "clearing old data"
    del asocL
#    del asocR
#if 'bL' in locals():
#    del bL
#    del bR
#    reload(ru)

asocL = ru.Brule(len(mod_names),len(col_names)-1,col_names[1:])
#bL = ru.Brule(len(mod_names),len(col_names)-1,col_names[1:])
asocL2 = ru.Brule(len(mod_names),len(col_names)-1,col_names[1:]) 
X=binary
    
asocL.trainAR(rPvalue,highT=1e-5,lowT=0.05)
#bL.trainFS(ratios[col_names[1:]])
asocR = asocL.predictAR(X,norm=False,sClass='flowerS')
#bR=bL.predictAR(X,norm=True,sClass='flowerS')

asocL2.trainAR2(rPvalue)
asocR2 = asocL2.predictAR(X,norm=False,sClass='flowerS')

temp = ["",'flowerS']
col_ind = [-1,-3]
#temp = col_names
# Evaluate results using prValue

y = []
yy =[]
classes = []
for i in range(1,len(temp)):
    glst = gLsts[i]
    name = temp[i]
    classes.append(name)
    yy.append(asocR[name].tolist())
    y.append([1 if x in glst else 0 for x in X.index])


vis.prCurve(y,yy,ratios,classes)
    


# End of PR curve

plt.figure()
for i in range(1,len(temp)):
    print i
    glst = gLsts[col_ind[i]]
    name = temp[i]
    scores =asocR[name]
    y= [1 if x in glst else 0 for x in scores.index]
    scores=scores.tolist()
    fpr,tpr,thresholds=metrics.roc_curve(y,scores)
    print fpr,tpr,thresholds
    plt.plot(fpr, tpr, label='My ROC curve for %s'%(name))
    
"""
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="upper left")
plt.show()
"""


for i in range(1,len(temp)):
    print i
    glst = gLsts[col_ind[i]]
    name = temp[i]
    scores =asocR2[name]
    y= [1 if x in glst else 0 for x in scores.index]
    scores=scores.tolist()
    fpr,tpr,thresholds=metrics.roc_curve(y,scores)
    print fpr,tpr,thresholds
    plt.plot(fpr, tpr, label='My ROC curve 2 for %s'%(name))



#plt.figure()
for i in range(1,len(temp)):
    print i
    glst = gLsts[col_ind[i]]
    name = temp[i]
    scores =bR[name]
    y= [1 if x in glst else 0 for x in scores.index]
    scores=scores.tolist()
    fpr,tpr,thresholds=metrics.roc_curve(y,scores)
    print fpr,tpr,thresholds
    plt.plot(fpr, tpr, label='ROC curve for %s'%(name))
    

plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="upper left")
plt.show()


