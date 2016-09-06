"""
This module is used to read the arabidopsis network 
and store it into igraph object

"""

import igraph as ig
import matplotlib.pyplot as plt

path = './at_net/reg_net.txt'


def getNet(path):

    name2id = dict()
    id2name = dict()

    edges = []

    ind = 0
    
    with open(path,'rb') as f:
        for line in f:
            tokens = line.rstrip().split('\t')
            if len(tokens)<13:
                continue
            source,target,istf,eType=tokens[0].upper(),tokens[3].upper(),tokens[5],\
                                      tokens[8]

            # Starting of the name <---> integet mapping
            if source in name2id:
                source = name2id[source]
            else:
                name2id[source]=ind
                id2name[ind]=source
                source = ind
                ind+=1

            if target in name2id:
                target = name2id[target]
            else:
                name2id[target]=ind
                id2name[ind]=target
                target =ind

                ind+=1

            # ===================================
            edges.append((source,target))

        f.close()

    return (ig.Graph(edges,directed=True),name2id,id2name)
            
                
def getDegreeHisto(lst,net,name2id):
    """
    This function returns the degree distribution of 
    the lists of genes (names)

    """

    indegree,outdegree,degree= net.indegree(),net.outdegree(),net.degree()

    indLst = [name2id[x] for x in lst if x in name2id]



    return ([indegree[i] for i in indLst],[outdegree[i] for i in indLst]\
            ,[degree[i] for i in indLst])



def plotDegreeHist(lst ,net,name2id):

    indegree,outdegree,degree = getDegreeHisto(lst, net,name2id)

    
    f,(ax1,ax2,ax3)=plt.subplots(3,sharex=True,sharey=True)
    
    
    ax1.hist(indegree)
    ax2.hist(outdegree)
    ax3.hist(degree)

    f.subplots_adjust(hspace=0)
    plt.show()
    

def getAverageDistance(lst,net,name2id):

    indLst = [name2id[x] for x in lst if x in name2id]

    vs = net.vs.select(indLst)

    return net.subgraph(vs).average_path_length()
    


    
