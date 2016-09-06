# This script should be used to test the ratio changes
# between different gene lists
# we can "execfile" run this script after loadFiles.py

import analysis as ana

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

print "Searching for frequent itemsets"

result = ana.compareFSS(dfLst,2)
lens = [x.shape[0] for x in dfLst]
ratios,rDiff,rPvalue,patterns = ana.createFssDF(result,mod_names,col_names,0,lens)

pPtns,nPtns = ana.findPatterns(rPvalue,patterns)







