# This module is used to visualize the histone modification
# enrichment level

import numpy as np
import matplotlib.pyplot as plt


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
        T= T[T[l]==1]

    for l in lno:
        T= T[T[l]==0]

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



def drawHistone(data,title,xs=range(-1000,1100,100),sf=None):
    """
    The input data is a dict of np arrays.
    Where the keys are the modification names
    and values is the location-specific enrichment levels
    """

    if all([x.shape[0]==0 for x in data.values()]):
           return 
    
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
    plt.xlim([-1200,1200])
    plt.ylim([-1,50])
    plt.xlabel('Location')
    plt.ylabel('Enrichment Level')
    plt.legend()
    plt.title(title)
    if sf is None:
        plt.show()
    else:
        plt.savefig(sf)
        



        
        
    

    
