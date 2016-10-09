# This module has the class for wrapping the
# histone modification data
import analysis as ana
import pandas as pd
import numpy as np
from scipy.stats import binom


class Mdata:
    """
    Mdata is the encapsule for the histone modificaiton
    data where the only binary data is stored for now.
    we also add the labels
    """


    
    def __init__(self,binary,mod_names):
        self.binary=binary
        self.mods= mod_names
        self.N = binary.shape[0]
        self.M = len(mod_names)
        self.labels = pd.DataFrame(index=binary.index)
        self.ratios = None
        self.rPvalue = None
        self.classes=[]
        
    def addLabels(self,labels,lName):
        """
        labels is a Series,list or dataframe
        of one column
        lName is the label name :flowering, defense
        """
        self.labels[lName] = 0
        self.labels.loc[labels,lName]=1
        self.classes.append(lName)
        
        
    def calcAllRatios(self):
        """
        After calling addLabels and findPatterns
        we calculate the ratios,rpvalue,and etc
        """
        patterns =[]
        self.lens = self.labels.sum()
        df = pd.DataFrame(index=self.pgMap.columns)
        for col in self.labels.columns:
            for ind in df.index:
                df.loc[ind,col]=len(\
            self.binary[(self.pgMap[ind]==1) & \
            (self.labels[col]==1)])/float(self.lens[col])
        self.ratios=df


        
    def calcrPvalue(self,base='all'):
        if self.ratios == None:
            self.calcAllRatios()
        df = pd.DataFrame(index=self.ratios.index)
        for col in self.labels.columns:
            print col
            temp = self.ratios.loc[:,col]\
                   *self.lens[col].astype(int)
            df[col] = binom.cdf(temp,\
                self.lens[col],self.ratios[base])
        self.rPvalue = df
    

        
    
    def findPatterns(self,min_support):
        """
        Find the patterns with the min_support
        """
        
        
        report = ana.xfindFS(self.binary,min_support)
        temp=ana.findFSGenes(report,self.binary)
        self.pgMap = pd.DataFrame(index=self.binary.index)

        for p in temp:
            self.pgMap[p]=0
            self.pgMap.loc[temp[p],p]=1

    def findPScores(self,mode='ar',ratios=None):
        """
        """
        if self.rPvalue==None:
            self.calcrPvalue()
            
        if mode == 'ar':
            if ratios is None:
                self.pScore= np.log10(self.rPvalue)-np.log10(1-self.rPvalue)
            else:
                self.pScore=np.log10(self.rPvalue).multiply(ratios['all'],axis='index')-np.log10(1-self.rPvalue)*ratios

    
    def predictALL(self,mode='norm'):
        """
        reevaluate all genes
        """
        result =None
        m = self.pgMap
        p = self.pScore
        if mode =='norm':
            result = self.binary.apply(lambda row:\
        p.ix[m.columns[m.ix[row.name].nonzero()]].sum() if len(m.ix[row.name].nonzero()[0])>0 else 0,axis=1)
        elif mode == 'maxmin':
            
            result = self.binary.apply(lambda row:\
        p.ix[m.columns[m.ix[row.name].nonzero()]].max()+p.ix[m.columns[m.ix[row.name].nonzero()]].min()if len(m.ix[row.name].nonzero()[0])>0 else 0,axis=1)

        return result
          
    
