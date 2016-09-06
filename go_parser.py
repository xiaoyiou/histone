"""A
This module deals with the GO terms and provide some
basic topic extraction analysis
"""

from yxtools import dictInsert
from nltk.tokenize import RegexpTokenizer
from stop_words import get_stop_words
from nltk.stem.porter import PorterStemmer

from collections import Counter

import matplotlib.pyplot as plt
import pandas as pd


def goReader(path):
    """
    This function creates two dictionaries
    1. rel between genes and GO terms ID as pandas DataFrame

    2. rel between GO term ID and descriptions a dictionary
    """

    gene =[]

    kind =[]

    relation = []

    goID = []

    goDesc = []

    slim=[]

    aspect=[]
    
    with open(path,'rb') as f:
        for line in f:
            tokens = line.rstrip().split('\t')

            gene.append(tokens[0])
            kind.append(tokens[1].split(':')[0])
            relation.append(tokens[3])
            goID.append(tokens[5])
            goDesc.append(tokens[4])
            aspect.append(tokens[7])
            slim.append(tokens[8])
            
        f.close()

    temp= pd.DataFrame({'gene':gene,'kind':kind,'rel':relation\
                         ,'goid':goID,'godesc':goDesc,'slim':slim,'aspect':aspect})

    return temp.set_index('gene')

def findGlst(gnames,mod_names):
    return frozenset([mod_names.index(x) for x in gnames])



def findDist(go,report,pattern,col,mod_names):
    
    glst =report[findGlst(pattern,mod_names)] 

    
    plt.figure()

    ax = go[go.index.isin(glst)][col].value_counts().plot(kind='barh')
    
#    for c in ax.containers:
#        plt.setp(c,width=1)

    plt.gcf().subplots_adjust(left=.35)
    
    plt.tight_layout()

    plt.title('totla number of GO terms: '+str(go[go.index.isin(glst)].shape[0]))

    plt.show()
    



    

            
"""
def findTopics(go,report,pattern,mod_names,top=50):

    glst =report[findGlst(pattern,mod_names)] 

    doc_set = go[go.index.isin(glst)]['godesc'].tolist()

    
    doc = ' '.join(doc_set)

    tokenizer = RegexpTokenizer(r'\w+')
    
    tokens = list(set(tokenizer.tokenize(doc)))

    en_stop = get_stop_words('en')

    stopped_tokens= [i for i in tokens if not i in en_stop]
    
    p_stemmer= PorterStemmer()
    
    texts = [p_stemmer.stem(i) for i in stopped_tokens]

    
    textx = unicodedata.normalize('NFKD', ' '.join(texts)).encode('ascii','ignore')

    

    counts =Counter(textx.split())

    lst = counts.most_common(top)

    for (x,y) in lst:
        print x,': ',y 
"""


def findGoTerms(go,pattern,report,mod_names):

    glst = report[findGlst(pattern, mod_names)]
    
    return go[(go.index.isin(glst)) & (go['rel']=='involved in')]
