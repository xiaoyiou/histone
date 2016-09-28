# This script is to test the expanding idea
import vis
import rule as ru
import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt
from functools import partial
import analysis as ana

def expandAR(mod_names,col_names,dfLst,glst,N=5,cname='flowerS',innerT =0, outterT =0,dfInd=-3):
    """
    N is the max number of iterations allowed
    binary : the binary input data
    labels : the 1/0 class labels
    thresh : the threshold for reliable genes
    """

    
    learner = ru.Brule(len(mod_names),len(col_names)-1,col_names[1:])


    X = dfLst[0]
    n = X.shape[0]
    fprs = []
    tprs = []
    nums = []

    for i in range(N):
        m = len(glst)
        print "Iteration {0:d} with {1:d}".format(i,m)
        dfLst[dfInd]=X.ix[glst]
        lens = [x.shape[0] for x in dfLst]
        print lens
        rr, _ = ana.compareFSS(dfLst,2,mod_names)
        ratios,rDiff,rPvalue=ana.createFssDF(rr,mod_names,col_names,0,lens)
        learner.trainAR(rPvalue,ratios=ratios)
        result = learner.predictAR(X,norm=True,sClass=cname)
        scores = result[cname]
        print "Average Inner Score", (np.mean(scores.ix[glst]))
        olst = scores[(scores>outterT) &\
                    (~scores.index.isin(glst))].index.tolist()
        print " {0:d}  are rated high"\
            .format(len(olst))

        nlst = scores[(scores>innerT) &\
                      (scores.index.isin(glst))].index.tolist()
        print " {0:d} From {1:d} labels are rated high"\
            .format(len(nlst),m)

        nlst = list(set(nlst).union(olst))
        nums.append(len(nlst))
        print "new label is of size {0:d}".format(nums[-1])
        if nlst == glst:
            print "same labels : stopping"
            break
        y= [1 if x in glst else 0 for x in scores.index]
        print sum(y)
        fpr,tpr,thresholds = metrics.roc_curve(y,scores)
        fprs.append(fpr)
        tprs.append(tpr)
        glst = nlst

    return (fprs,tprs,nums,learner)
        


def expandFS(mod_names,col_names,dfLst,glst,N=5,cname='flowerS',innerT =0, outterT =0,dfInd=-3):
    """
    N is the max number of iterations allowed
    binary : the binary input data
    labels : the 1/0 class labels
    thresh : the threshold for reliable genes
    """



    X = dfLst[0]
    n = X.shape[0]
    fprs = []
    tprs = []
    nums = []

    for i in range(N):
        
        learner = ru.Brule(len(mod_names),len(col_names)-1,col_names[1:])
        m = len(glst)
        print "Iteration {0:d} with {1:d}".format(i,m)
        dfLst[dfInd]=X.ix[glst]
        lens = [x.shape[0] for x in dfLst]
        print lens
        rr, _ = ana.compareFSS(dfLst,2,mod_names)
        ratios,rDiff,rPvalue=ana.createFssDF(rr,mod_names,col_names,0,lens)
        learner.trainFS(ratios[col_names[1:]])
        result = learner.predictAR(X,norm=True,sClass=cname)
        scores = result[cname]
        print "Average Inner Score {0:.3f}".format(np.mean(scores.ix[glst]))
        olst = scores[(scores>outterT) &\
                    (~scores.index.isin(glst))].index.tolist()
        print " {0:d} are rated high".format(len(olst))

        nlst = scores[(scores>innerT) &\
                      (scores.index.isin(glst))].index.tolist()
        print " {0:d} From {1:d} labels are rated high"\
            .format(len(nlst),m)

        nlst = list(set(nlst).union(olst))
        nums.append(len(nlst))
        print "new label is of size {0:d}".format(nums[-1])
        if nlst == glst:
            print "same labels : stopping"
            break
        y= [1 if x in glst else 0 for x in scores.index]
        fpr,tpr,thresholds = metrics.roc_curve(y,scores)
        fprs.append(fpr)
        tprs.append(tpr)
        glst = nlst
        del learner
    return (fprs,tprs,nums,learner)
