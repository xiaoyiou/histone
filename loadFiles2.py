"""
This script loads the modifications and genome file for 
arabidopsis in the seedling stage

"""

import os,sys
import cPickle as P
import os.path
import analysis as ana
import pandas as pd
import numpy as np

aPath='T8_TSS1k.bed'
#bPath='baseline2.p'

prefix = '10day/bgs/'


tPath='temp2.csv'
ttPath='temp2'

chrlens={'chr1':30432562,
         'chr2':19705358,
         'chr3':23470804,
         'chr4':18585041,
         'chr5':26992727,
         'chrC':154477,
         'chrM':366923}


mod_files = ['5mC.bed','H3.bed','H3K27me3.bed','H3K4me2.bed','H2BUb.bed','H3K27me1.bed'\
,'H3K36me3.bed','H3K4me3.bed']

mod_names=['5MC','H3','H3K27ME3','H3K4ME2','H2BUB','H3K27ME1','H3K36ME3','H3K4ME3']




if not os.path.isfile(tPath):
    sep=' '
    
    templst = [prefix+f for f in mod_files]

    command = 'bedtools intersect -wa -wb -loj -a ' +aPath+ ' -b ' +sep.join(templst)+ ' >' + ttPath
    
    os.system(command)

#    open(tPath,'w').writelines([ line for line in open(ttPath) if  line.split('\t')[4]!='.'])
     

############################################################################



data=ana.getMods(ttPath,countType=np.float)
gData = ana.getEnrich2(data,ana.getData(data))
threshs = ana.calcThreshAll(data,power=8)
binary =ana.maxBin(data,threshs)
