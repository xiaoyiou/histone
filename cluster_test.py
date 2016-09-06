"""
This script is used following loadFiles to conduct seriese of 
grouping schema for genes. 

1. Frequent itemset: genes with same pattern are grouped together

2. 

"""

import analysis as ana
import go_parser as go
import atNet2 as at



print "Loading Go Term Files"

goterms = go.goReader('./ATH_GO_GOSLIM.txt')

print "Searching for frequent itemsets"

ana.findFS(binary, mod_names, min_support=10, maxReport=50)

print "Getting report for gene lists"

report = ana.findFSGenes(ana.__findFS(binary,10),mod_names,binary)


print "Loading Gene Regulatory Network"

G,name2id,id2name = at.getNet(at.path)


pos = at.gt.sfdp_layout(G)

