# This script is to test
import vis
import rule as ru
import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt
if 'learner' in locals():
    del learner
    reload(ru)

learner = ru.Brule(len(mod_names),len(col_names)-1,col_names[1:])

X =binary

learner.trainAR(rPvalue,lowT=0.05,highT=0.05)
result = learner.predictAR(X,sClass='flowerS',norm=True)
#result = learner.predictAR(X,norm=True)

temp = ["",'flowerS']
#temp = col_names
# Evaluate results using prValue
"""
y = []
yy =[]
classes = []
ratios =[]
for i in range(1,len(temp)):
    glst = gLsts[i]
    name = temp[i]
    classes.append(name)
    yy.append(result[name].tolist())
    y.append([1 if x in glst else 0 for x in X.index])


vis.prCurve(y,yy,ratios,classes)
    


# End of PR curve





"""
plt.clf()
for i in range(1,len(temp)):
    print i
    glst = gLsts[i]
    name = temp[i]
    scores = result[name]
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


