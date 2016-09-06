# This script is only used to get the housekeeping gene list

path ='./housekeepers'


i=0
with open(path,'rb') as f:
    for line in f:
        tokens = line.rstrip().split('\t')
        gene = tokens[1]
        print gene.upper()
        i+=1
    f.close()

print i
    
