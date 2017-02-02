# This module is used to visualize the histone modification
# enrichment level

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd




def prepareData(gdata,mods,mod_names,lyes,lno,ldata,bNum=21):
    """
    labels is set of labels 
    and mods will be the set of mod_names
    It will return a dictionary where keys are the mod names
    and value is the numpy.array of actual enrichment levels
    """
    result = dict()
    
    # prepare labels:
    L = [mod_names.index(x) for x in mods]

    T = ldata.copy() # mdata.labels

    for l in lyes:
        T= T[T[l] == 1]

    for l in lno:
        T= T[T[l] == 0]

    glst = T.index.tolist()

    for i in L:
        modvs= []
        for g in glst:
            if len(gdata[(g,i)])!=bNum:
                continue
            modvs.append(gdata[(g,i)])
        result[mod_names[i]]=np.array(modvs,dtype=np.float)
        
    print "number of genes: " , len(glst)
    print "mod names:" ,mods
    return result



def drawHistone(data,title,xs=range(-1000,1100,100),sf=None,append = False):
    """
    The input data is a dict of np arrays.
    Where the keys are the modification names
    and values is the location-specific enrichment levels
    """

    if all([x.shape[0]==0 for x in data.values()]):
           return 
    if not append:
        plt.figure(figsize=(10,10))
    for mod in data:
        d = data[mod]

            
        dmean= d.mean(0)
        dstd= d.std(0)
        dupper=dmean+dstd
        dlower=np.maximum(dmean-dstd,1)
        print dmean,dstd,dupper,mod
        mplot=plt.plot(xs,dmean,label=mod)
        plt.fill_between(xs,dlower,dupper,alpha=0.3,\
                         color=mplot[0].get_color())

    if not append:
        plt.xlim([-1200,1200])
        plt.ylim([-1,50])
        plt.xlabel('Location')
        plt.ylabel('Enrichment Level')
        plt.title(title)
    plt.legend()
    
        
    if sf is None:
        plt.show()
    else:
        plt.savefig(sf)
        

def prepareDataLst(gdata,mods,mod_names,glst,bNum=21):
    """
    Use genelist create the subset of data
    """
    result = dict()
    # prepare labels:
    L = [mod_names.index(x) for x in mods]

    for i in L:
        modvs= []
        for g in glst:
            if len(gdata[(g,i)])!=bNum:
                continue
            modvs.append(gdata[(g,i)])
        result[mod_names[i]]=np.array(modvs,dtype=np.float)
        
    print "number of genes: " , len(glst)
    print "mod names:" ,mods
    return result

        

        
        
    
def drawHistoneLst(data, title, mod, xs=range(-1000,1100,100), sf = None, append=False, smooth=True, ymax = 50):
    """
    Here, we draw enrichment level for one modifcation
    across all genes in the data. No average, No std.
    """
    if all([x.shape[0]==0 for x in data.values()]):
           return 
    if not append or plt.fignum_exists(1):
        plt.figure(num=1,figsize=(10,10))
    
    d = data[mod]
    for row in d:
        if smooth:
            wd = 1
            for i in range(len(row)):
                row[i]=np.mean(row[max(0,i-wd):\
                min(len(row),i+wd)])
        mplot=plt.plot(xs,row,label=mod)

    if not append:
        plt.xlim([-1200,1200])
        plt.ylim([-1,ymax])
        plt.xlabel('Location')
        plt.ylabel('Enrichment Level')
        plt.title(title)
 
    
        
    if sf is None:
        plt.show()
    else:
        plt.savefig(sf)


def drawGenes(data, title, mod, xs=range(-1000,1100,100), sf = None, append=False, smooth=True, ymax = 50, fill=True):
    """
    Here, we draw different genes from different categories.
    data is a dictionary of np arrays
    {category: dict of mods}
    """
    mean_data ={}
    if all([len(x)==0 for x in data.values()]):
           return
    if not append or plt.fignum_exists(1):
        plt.figure(figsize=(10,10))

    for label in data:
        d = data[label][mod]
        dmean= d.mean(0)
        mean_data[label] = dmean
        dstd= d.std(0)
        dupper=dmean+dstd
        dlower=np.maximum(dmean-dstd,1)
        mplot=plt.plot(xs,dmean,label=label)
        if fill:
            plt.fill_between(xs,dlower,dupper,alpha=0.3,\
                         color=mplot[0].get_color())

    if not append:
        plt.xlim([-1200,1200])
        plt.ylim([-1,ymax])
        plt.xlabel('Location')
        plt.ylabel('Enrichment Level')
        plt.title(title)
        plt.legend()


    if sf is None:
        plt.show()
    else:
        plt.savefig(sf)

    return pd.DataFrame(mean_data)



    
    
    
