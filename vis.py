# This module is used to visualize the results from rule learning
# algorithms
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve as pr
from sklearn.metrics import average_precision_score as apr
from sklearn import metrics
from ggplot import *
import seaborn as sns
import pandas as pd
def prCurve(y_test,y_score,ratios,classes):
    """
    ;param y_test 
    ;param y_score 
    ;param ratio 
    ;param classes [] is the list of classes names
    """

    N = len(y_test)

    plt.figure()
    
    for i in range(N):
        y = y_test[i]
        yy = y_score[i]
        name = classes[i]
        prec, recall,_ = pr(y,yy)

        avpr = apr(y,yy)

        plt.plot(recall,prec,label='%s: APR=%.2f'%(name,avpr))
        

    plt.xlim([0,1])
    plt.ylim([0,1.05])

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title("PR curve for multi-class")
    plt.legend(loc='top right')
    plt.show()
        

def roc(fprs,tprs,names):
    """
    Plot roc curves with different names
    using ggplot
    """

    for i in range(len(names)):
        auc = metrics.auc(fprs[i],tprs[i])
        plt.plot(fprs[i], tprs[i], label=\
                 '%dth iter %s with AUC=%.2f'%(i,names[i],auc))
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="upper left")
    plt.show()



def visPtns(scores,binary,glst,target,step=100):
    temp = scores.ix[glst].sort_values(by=target,ascending=False).index

    data = pd.DataFrame(columns=binary.columns)
    i = 0
    N = temp.shape[0]
    nData = binary.ix[temp]
    
    while i<N:

        x = nData.iloc[range(i,min(i+step,N))].mean()
        
        data.ix[str(i)]=x
        i+=step
    plt.clf()
    sns.heatmap(data)

    
    
    
    

    
