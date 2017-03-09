# this module should take care of all weird stuff
from yxtools import getTabData
import sys

def getPatternInfo(data, path, obj, col_inds=[1, 2, 5, 6, 7], k=10, opath=None):
    """

    Parameters
    ----------
    data: should be the sorted saxGrps
    path: should be the path of gene names
    obj: the object of cdata.MData

    Returns: a list of patterns with sorted genes based on function labels
    -------

    """
    stdout = sys.stdout
    if opath:
        sys.stdout = open(opath, 'wb')
    names = dict()
    labels = obj.labels
    cols = labels.columns.tolist()
    col_names = [obj.classes[x] for x in col_inds]
    names = getTabData(path, delim='\t')
    res = []
    for pattern, glst in data:
        if len(glst) < k:
            break
        print '#'*50
        print pattern
        print '#' * 50 + '\n'
        t = labels.ix[glst]
        grps = t.groupby(by=col_names).groups
        for group, lst in grps.items():
            name = []
            for i,v in enumerate(group):
                if v != 0:
                    name.append(col_names[i])
            print '$' * 50
            print 'Labels: %s'%(' '.join(name) if len(name) else 'NO')
            for gene in lst:
                print gene, names.get(gene)
            print '$' * 50
            print ""
        print ''
    if opath:
        sys.stdout.close()
        sys.stdout = stdout
