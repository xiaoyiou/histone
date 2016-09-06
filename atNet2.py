"""
This module is the same as the atNet module
but we use graph_tool instead of igraph here

"""

import graph_tool.all as gt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde 
from sklearn import datasets, linear_model


path='./at_net/reg_net.txt'

def getNet(path, dir=True):
    name2id = dict()
    id2name = dict()
    edges = []
    ind = 0
    
    with open(path,'rb') as f:
        for line in f:
            tokens = line.rstrip().split('\t')
            if len(tokens)<13:
                continue
            source,target,istf,eType=tokens[1].upper(),tokens[4].upper(),\
                                      tokens[5],tokens[8]

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

    g = gt.Graph(directed=dir)
    g.add_edge_list(edges)
    return (g,name2id,id2name)


def drawEntireNet(g,pos=None, glst=None,nmap=None,hc=[1,0,0,1]):
    """
    drawEntireNet draws the entire network with highlighted
    vertexes.
    """
    lpos=None
    if pos!=None:
        lpos=pos
    else:
        lpos=gt.sfdp_layout(g)
        
    state = g.new_vertex_property('vector<double>')
    size = g.new_vertex_property('int')
    deg = g.degree_property_map('out')
    deg.a=4*(np.sqrt(deg.a)*0.5+0.4)
    
    for v in g.vertices():
        state[v]=[1,1,1,1]
        size[v]=1
        
    if glst!=None:
        idLst=[nmap[x] for x in glst if x in nmap]
        for ind in idLst:
            state[ind]=hc
            size[ind]=10
    
    gt.graph_draw(g, pos=lpos, vertex_fill_color=state,\
                  edge_color=[0.6,0.6,0.6,1],vertex_size=size)

def createDegreeMods(binary,id2name,G,nMods):
    """
    createDegreeMods creates the panda data frame 
    for further analysis 
    """

    dgrAll = G.degree_property_map('total')

    dgrIn = G.degree_property_map('in')

    dgrOut = G.degree_property_map('out')

    pdNames = []

    N = len(dgrIn.a)

    for i in range(N):
        if i in id2name:
            pdNames.append(id2name[i])
        else:
            pdNames.append('NULL')

    degs = pd.DataFrame({'id':range(N),'dgrA':dgrAll.a.tolist()\
                         ,'dgrI':dgrIn.a.tolist(),'dgrO':dgrOut.a.tolist(),\
                         'name':pdNames})
            
    tmp = binary.merge(degs,left_index=True,right_on='name',how='inner')

    tmp['modN']=tmp[range(nMods)].sum(axis=1).tolist()

    return tmp

def visDegModRel(df, nMods, lst =None):
    """
    findDeggModRel try to find the relationship between # of mods
    and degree of genes in a global sense, if lst is set to None
    otherwise we only consider the genes in the lst (intersection)
    """

    texts = ['Total', 'In', 'Out']
    
    ys = [df[x].tolist() for x in ['dgrA','dgrI','dgrO']]

    x = df['modN']

    
    
    ii= 0

    fig =plt.figure()
    
    
    for y in ys:
        plt.figure()
        #H,xedges,yedges = np.histogram2d(x,y)
        #X,Y =
        plt.hist2d(x,y)
        plt.colorbar()
        plt.show()
        ii+=1
        
    

def findDegModRel(data, nMods,norm=False,filt=0,lst =None):
    """
    findDegModRel function tries to find the linear correlation between
    the degree and number of modifications
    """

    if filt>0:
        data = data[data['dgrA']>filt]
    

        
        
    regr = linear_model.LinearRegression()

    y1,x1 = data['dgrA'],data['modN']

    if norm:
        y1 = (y1-y1.mean())/(y1.max()-y1.min())
        x1 = (x1-x1.mean())/(x1.max()-x1.min())
    
    x,y = np.array([[i] for i in x1]),np.array([[j] for j in y1])
    
    xy = np.vstack([x1,y1])
    z = gaussian_kde(xy)(xy)    
    regr.fit(x,y)

    print 'Coefficient',regr.coef_
    print 'Residual sum of squares (training): %f'\
        %np.mean((regr.predict(x)-y)**2)

    plt.plot(x,regr.predict(x),color='red',linewidth=3)

    plt.scatter(x,y,c=z,s=250,edgecolor='')


    plt.show()


    
    
def findCores(G,glst,nmap,thresh=3, directed=False):
    """
    findCores will create the subnetwork 
    """

    idMap = dict()
    g = gt.Graph(directed=directed)
    idLst=[nmap[x] for x in glst if x in nmap]

    for ind in idLst:
        v = g.add_vertex()
        idMap[v] = ind

    dist=gt.graph_tool.topology.shortest_distance(\
            G,directed=directed)
    
    filt = []    
    for x in g.vertices():
        source = idMap[x]
        temp = dist[source]
        for y in g.vertices():
            target = idMap[y]
            if (not directed) and (target,source) in filt:
                continue
            else:
                filt.append((source,target))
            if source ==target:
                continue
            tmp = temp[target]
            if tmp > 1e8:
                continue
            else:
                if tmp<thresh:
                    g.add_edge(x,y)
                else:
                    continue

    pos = gt.sfdp_layout(g)
    gt.graph_draw(g,pos=pos)

    
def relatePattern(report, name2id):
    """
    relatePattern creates the dictionary for
    id -> pattern mapping. we only keep the modid 
    """
    result = dict()
    for key in report:
        lst = report[key]
        for gene in lst:
            if gene in name2id:
                result[name2id[gene]]=key
    return result

def __DFS(G,source,pattern,color):
    result = []
    edges = []
    color[source]=1
    v = G.vertex(source)
    result.append(source)
    for i in v.all_neighbours():
        key = int(i)
        if not key in pattern or not source in pattern:
            continue
        if color[key]==0 and pattern[key].issubset(pattern[source]):
            edges.append((source,int(i)))
            te,tr=__DFS(G,key,pattern,color)
            edges+=te
            result+=tr
    return edges,result

def travModsPaths(G, source, pattern):
    """
    pattern is the id->pattern dictionary from relatePattern
    This method does a DFS search with contraint of 
    modifications 
    """
    
    color = [0]*G.num_vertices()    
    edges, result = __DFS(G,source,pattern,color)
    edges,result =list(set(edges)),list(set(result))
                
    return edges,result


def highLightSubNet(G,gIDlst,pos=None,hc=[1,0,0,1],size=10):

    lpos=None
    if pos!=None:
        lpos=pos
    else:
        lpos=gt.sfdp_layout(G)

    state = G.new_vertex_property('vector<double>')

    for v in gIDlst:
        state[v]=hc
    
    gt.graph_draw(G,pos=lpos,vertex_fill_color=state,size=size)


def highLightSubNets(G,glst,nmap,pattern,pos=None):

    setGes=dict()
    
    lpos=None
    if pos!=None:
        lpos=pos
    else:
        lpos=gt.sfdp_layout(G)

#    state = G.new_vertex_property('vector<double>')
    state = G.new_vertex_property('int')
    size  = G.new_vertex_property('int')
    prop = G.new_edge_property('int')
    width = G.new_edge_property('int')
    
    for v in G.vertices():
        state[v]=0
        size[v]=.1
        
    for e in G.edges():
        prop[e]=0
        width[e]=.1
        
    ii=0
    
    for g in glst:
        if not g in nmap:
            continue
        gId = nmap[g]
        edges,reachable = travModsPaths(G,gId,pattern)
        
        for i in reachable:
            state[i]=10
            size[i]=10
            
        for i in edges:
            if i in G.edges():
                prop[i]=ii
                width[i]=10
        ii+=1
        
    gt.graph_draw(G,pos=lpos,vertex_fill_color=state,vertex_size=size,\
                  edge_color=prop,edge_pen_width=width)

    
