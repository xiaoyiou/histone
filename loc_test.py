# This is a sanity check of the hypothesis for importance
# of position and raw enrichment level

import vish as vh
import numpy as np

stop_mods = ['H3','H3K9ME2','H3K27ME1']
stop_labels=['flowerS','flowerB']
resultPath='/Users/xiaoyiou/OneDrive/result/Nov_13/all/'
mods = mod_names
#mods= [x for x in mod_names if not x in stop_mods]
#mods=['H3K4ME3']
#labels=[x for x in young.classes[1:] if not x in stop_labels]
labels=['flowerN', 'stress']

"""
for x in labels:
    for y in labels:
        if x==y:
            continue
        title = x+'-'+y
        d = vh.prepareData(gData,mods,mod_names,\
                           [x],[y],young.labels)
#        vh.drawHistone(d,title,sf='/Users/xiaoyiou/OneDrive/result/Nov_1/'+title+'.png')
        vh.drawHistone(d,title)
"""
cor_data = {'pearson':[],'spearman':[],'norm':[],'euclidean':[]}
cor_indxs = []
all_data = {}
for l in labels:

    all_data[l] = vh.prepareDataLst(gData,mods,mod_names,geneLsts[l])
for m in mods:
    title = "Distribution of modification: %s" % m
    mean_data = vh.drawGenes(all_data,title,m,ymax=50, fill=True,sf=resultPath+m+'.png')
    cor_indxs.append(m)
    cor_data['pearson'].append(mean_data['stress'].corr(mean_data['flowerN'],method='pearson'))
    cor_data['spearman'].append(mean_data['stress'].corr(mean_data['flowerN'],method='spearman'))
    cor_data['norm'].append(np.linalg.norm(mean_data['stress']-mean_data['flowerN'],ord =1))
    cor_data['euclidean'].append(np.linalg.norm(mean_data['stress']-mean_data['flowerN'],ord =2))

pd.DataFrame(cor_data,index=cor_indxs).to_csv(resultPath+'all.csv')
