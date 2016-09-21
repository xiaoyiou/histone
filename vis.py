# This module is used to visualize the results from rule learning
# algorithms
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve as pr
from sklearn.metrics import average_precision_score as apr

def prCurve(y_test,y_score,ratios,classes):
    """
    ;param y_test 
    ;param y_score 
    ;param ratio 
    ;param classes [] is the list of classes names
    """

    N = len(y_test)

    plt.clf()
    
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
        
        