# This script should be used to test the ratio changes
# between different gene lists
# we can "execfile" run this script after loadFiles.py

import analysis as ana
import rule as ru
from yxtools import dictInsert





glst_paths = ['defense.glst','development.glst'
              ,'flower.glst','flowering.glst','stress.glst'\
              ,'stimulus.glst']

inds = [0,1,0,0,1,1]

col_names = ['all','defense','develop','flowerB'\
             ,'flowerS','stress','stimulus']
prefix ='genelst/'

dfLst = [binary]+[ana.selectDataGL\
         (binary,prefix+glst_paths[i],inds[i])\
         for i in range(len(glst_paths))]


gLsts = [x.index.tolist() for x in dfLst]
print "Searching for frequent itemsets"

result,pGlst = ana.compareFSS(dfLst,2,mod_names) # pGlst is the dict to transform patterns to gene lists
lens = [x.shape[0] for x in dfLst]
ratios,rDiff,rPvalue = ana.createFssDF(result,mod_names,col_names,0,lens)


print "Generating postive and negative patterns"
pPtns,nPtns = ana.findPatterns(rPvalue,lowT=5e-2,highT=5e-2,ratios=ratios,threshP=0.05,threshN=0.05)



"""
print "start testing different thresholds"
performance = dict()
lowT=5e-2
highT=5e-2
M=1

for i in range(M):
    print i
    pPtns,nPtns = ana.findPatterns(rPvalue,lowT=lowT/(10**i),highT=highT/(10**i),ratios=ratios,threshP=0,threshN=0)

    print "start learning process"
    learner = ru.Brule(len(mod_names),len(col_names)-1)

    learner.maxPtnTrain(pPtns,nPtns)

    #learner.naiveTrain(pPtns,nPtns)
    
    labels =learner.predict(binary)

    perf = ru.evaluate(labels,gLsts,col_names)
   
    for key in perf:
        dictInsert(performance,key,perf[key])

print performance
"""
