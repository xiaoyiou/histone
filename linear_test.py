from sklearn import linear_model
from itertools import cycle
#from sklearn.model_selection import 


fig = plt.figure(facecolor='white',edgecolor='none')
linecycler= cycle(["-","--","-.",":","-v"])
m = young.pgMap.copy()

for name in young.classes[1:]:
    if name in ['flowerS','flowerB']:
        continue
    logreg = linear_model.LogisticRegression(C=1e5)
#    XX=X
#    XX = m.apply(lambda row: row*young.pScore[name],axis=1)

    XX = m.apply(lambda row: row*young.ratios[name],axis=1)
 
#    temp =X[young.labels[name]==1]
#    temp= temp.sum()/X.shape[0]
#    XX=X.apply(lambda row: row*temp,axis=1)

    
    y = young.labels[name]
    logreg.fit(XX,y)
    scores=logreg.predict_proba(XX)[:,1]
    fpr,tpr,thresholds=metrics.roc_curve(y,scores)
    plt.plot(fpr, tpr,next(linecycler), label='%s, AUC: %.2f'%(name,metrics.auc(fpr,tpr)))
    

plt.plot([0, 1], [0, 1], 'k--',label='Random classifier')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC for function specific binary classifier')
plt.legend(loc="lower right",fontsize=18)
plt.show()

fig.savefig('~/OneDrive/result/Oct_31/ratios.png',faceolor='white',edgecolor='none',transparent=True)
