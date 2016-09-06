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

aPath='T8_TSS2k.bed'
bPath='baseline2.p'

prefix = '10day/bgs/'

if len(sys.argv)>1:
    # The gene annotation bed file path
    aPath =sys.argv[1]
    
    # The baseline count for each modification
    bPath = sys.argv[2]
    
    # The temp file path

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

    open(tPath,'w').writelines([ line for line in open(ttPath) if  line.split('\t')[4]!='.'])
     
# load the files for calculating the average mod number in each bin
# result is stored in mod:chromosome dictionary
baseline = dict()
bdf=[]
if os.path.isfile(bPath):
    baseline = P.load(open(bPath,'rb'))
else:
    for i in range(len(mod_files)):
        fname=mod_files[i]
        with open(prefix+fname,'rb') as f:
            
            for line in f:
                tokens=line.rstrip().split()
                count=float(tokens[3])
                key=(i,tokens[0][:-1]+tokens[0][-1].upper())
                if key in baseline:
                    baseline[key]+=count
                else:
                    baseline[key]=count
            f.close()
    for (mod,chrom) in baseline:
        bdf.append((mod,chrom,baseline[(mod,chrom)]/float(chrlens[chrom]/100)))
    baseline=pd.DataFrame(bdf,columns=['modid','chrom','mean'])
    P.dump(baseline,open(bPath,'wb'))
    

############################################################################




data=ana.getMods(tPath,countType=np.float32)
maxes=ana.findMax(data)

binary=pd.merge(maxes,baseline)
binary['appear']=0
binary.ix[binary['count']>binary['mean']*3,'appear']=1
binary=binary.pivot(index='gene',columns='modid',values='appear')
