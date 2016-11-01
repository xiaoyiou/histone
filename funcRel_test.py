ps = young.pScore
pgMap = young.pgMap

max_ptns = dict()
labels = dict()
for i  in ps.columns[1:]:
    max_ptns[i]=ps[i].argmax()
    labels[i]=young.labels.index[young.labels[i].nonzero()].tolist()
# Now for each pattern calculate the precision.
"""
precs = dict()

for i in max_ptns:
    pattern = max_ptns[i]
    A = pgMap[[pattern]]
    lst = A[A.iloc[:,0]==1].index
    precs[i]=dict()
    for j in labels:
#        precs[i][j]=len(lst.intersection(labels[j]))/float(len(lst))
        precs[i][j]=len(lst.intersection(labels[j]))/float(len(labels[j]))
"""

for i in max_ptns:
    p1 = max_ptns[i]
    precs[i]=dict()
    for j in labels:
        p2= max_ptns[j]
        A = pgMap[[p1.union(p2)]]
        lst = A[A.iloc[:,0]==1].index
#        precs[i][j]=len(lst.intersection(labels[j]))/float(len(lst))
        precs[i][j]=len(lst.intersection(labels[j]))/float(len(labels[j]))


precs_df=pd.DataFrame.from_dict(precs,orient='index')
precs_df=precs_df.ix[precs_df.columns]

    
