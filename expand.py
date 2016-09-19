# This script is to test the expanding idea


print "start intial learning process"
learner = ru.Brule(len(mod_names),len(col_names)-1)

learner.maxPtnTrain(pPtns,nPtns,r=False,a=False)

labels =learner.predict(binary)


col2=['all','flowerS']
for i in range(10):
    ll = labels['flowerS']
    dfLst2 = [binary,binary.ix[ll[ll==1].index]]

    lens2 = [x.shape[0] for x in dfLst2]
    result2,pGlst2 = ana.compareFSS(dfLst2,2,mod_names) # pGlst is the dict to transform patterns to gene lists
    ratios2,rDiff2,rPvalue2 = ana.createFssDF(result2,mod_names,col2,0,lens2)


    print "Generating postive and negative patterns"
    pPtns2,nPtns2 = ana.findPatterns(rPvalue2,lowT=5e-2,highT=5e-2,ratios=ratios,threshP=0.05,threshN=0.05)


    print "start learning process"
    learner2 = ru.Brule(len(mod_names),len(col2)-1)

    learner2.maxPtnTrain(pPtns2,nPtns2,r=False,a=False)

    labels =learner2.predict(binary)

    for x in labels:
        print x, labels[x].sum()
