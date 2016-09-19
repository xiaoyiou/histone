# This modiule deals with the rule learning and testing
import numpy as np

def evaluate(labels,glsts,col_names):

    result = dict()



    for i in range(len(col_names)-1):
        key = col_names[i+1]
        t = glsts[i]
        p = labels[key][labels[key]==1].index
        tp = p.intersection(t)
        tpn= float(len(tp))
        if len(p)==0:
            p=[1]
        result[key] =(tpn/(len(p)),tpn/(len(t)))
        
    return result

        
class Brule(object):
    """
    Brule is the class for binary rules
    """
    classes = []
    
    pB = dict()
    # This is list of numpy arrays of 1/0s
    # to represent the postively related patterns
    nB = dict()
    # to represent the negatively related patterns
    
    def __init__(self,N,M):
        """
        N is the number of modifications
        M is the number of classes
        """
        self.N=N
        self.M=M

    def predictHelper(self,row,pb,nb):

        if any([all(row | ~x) for x in pb])\
        and not any([all(row | ~x) for x in nb]):
            return 1
        else:
            return 0

        
    def predict(self,df):
        #if not len(row)==self.N:
        #    raise AssertionError('pattern/data lengths are not equal')

        data = dict()
        df = df.astype('bool')
        for key in self.pB:            
            pb, nb=self.pB[key],self.nB[key]
            data[key]=df.apply((lambda x: \
                self.predictHelper(x,pb,nb)),axis=1)
            
        return data
    
    def __addRules(self,pPtn):
        """
        This module does nothing but
        adding the rules as binary reps to the learner
        """
        lst = []
        for ptn in pPtn:
            temp=[0]*self.N
            for i in ptn:
                temp[i]=1
            lst.append(np.array(temp,dtype='bool'))
        return lst
        

    def __naiveTrain(self,pPtn,nPtn,name):
        """
        This is a helper method for naiveTrain
        to create rules for a single class
        The algorithm is fairly simple.
        For all positive rules, we take the union of 
        postive mods and subtract the union of 
        negatively related mods
        """
        self.pB[name]=self.__addRules(pPtn)
        self.nB[name]=self.__addRules(nPtn)
        

    def findMaxPtns(self,ptns):
        """
        Here patterns are frozensets 
        we eliminate all the patterns with at least one
        superset
        """
        result = []
        for ptn in ptns:
            if not any([x.issuperset(ptn)\
                        for x in ptns if x!=ptn]):
                result.append(ptn)
        return result
    
        
    def naiveTrain(self,pPtns,nPtns):
        """
        The naive trainer will train the 
        """
        for key in pPtns:
            self.__naiveTrain(pPtns[key],nPtns[key],key)

            

    def maxPtnTrain(self,pPtns,nPtns,r=False,a=False):
        """
        (r)educe: is the flag for reducing conflicting
        (a)nnihlate: is the flag for removing identical p/n ptns
        """

        universe = frozenset(range(self.N))
        
        for key in pPtns:
            pb = self.findMaxPtns(pPtns[key])
            nb = self.findMaxPtns(nPtns[key])

            if r:
                for x in pb:
                    for y in nb:
                        if x==y:
                            if a:
                                pb.remove(x)
                                nb.remove(y)
                            else:
                                continue
                        if y.issuperset(x):
                            pb.remove(x)
                        if x.issuperset(y):
                            nb.remove(y)

            
                
            self.__naiveTrain(pb,nb,key)

        
