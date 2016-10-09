"""
This script loads the modifications and genome file for 
arabidopsis th

"""

import os,sys
import cPickle as P
import os.path
import analysis as ana
import pandas as pd



aPath='T8_TSS1k.bed'
bPath='baseline.p'

prefix = 'young/'

# The following part is commented out because notebook
# is passing arguments to the script. uncomment if
# not running in interactive notebook
"""
if len(sys.argv)>1:
    # The gene annotation bed file path
    aPath =sys.argv[1]
    
    # The baseline count for each modification
    bPath = sys.argv[2]
    
    # The temp file path
"""
tPath='temp.csv'


chrlens={'chr1':30432562,
         'chr2':19705358,
         'chr3':23470804,
         'chr4':18585041,
         'chr5':26992727,
         'chrC':154477,
         'chrM':366923}


mod_files=['gsm701927_h3k18ac.bedgraph','gsm701926_h3k9me2.bedgraph',\
'gsm701924_h3k4me3.bedgraph','gsm701930_h3k36me2.bedgraph',\
           'gsm701925_h3k9ac.bedgraph','gsm701923_h3k4me2.bedgraph',\
'gsm701931_h3k36me3.bedgraph','gsm701932_h3.bedgraph',\
'gsm701929_h3k27me3.bedgraph','gsm701928_h3k27me1.bedgraph']

mod_names=['H3K18AC','H3K9ME2','H3K4ME3','H3K36ME2','H3K9AC',\
           'H3K4ME2','H3K36ME3','H3','H3K27ME3','H3K27ME1']




if not os.path.isfile(tPath):
    sep=' '

    templst=[prefix+f for f in mod_files]

    command = 'bedtools intersect -wa -wb -loj -a ' +aPath+ ' -b ' +sep.join(templst)+ ' >' + tPath
    
    

    os.system(command)


     
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
                count=int(tokens[3])
                key=(i,tokens[0])
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

data=ana.getMods(tPath)
gData = ana.getEnrich2(data,ana.getData(data))
threshs = ana.calcThreshAll(data)
binary =ana.maxBin(data,threshs)

