# this module aims at creating GO term networks with frequent co-occurrence

import go_parser
import graph_tool.all as gt
import numpy as np
evidence = {'EXP':4, 'IDA':4,'IPI':4, 'IMP':4, 'IGI':4 ,'IEP':4,\
            'ISS':3, 'ISO':3, 'ISA':3, 'ISM':3, 'IGC':3, 'IBA':3,\
            'IBD':3, 'IKR':3, 'IRD':3, 'RCA':3, 'TAS': 2, 'NAS':2,\
            'IC':1, 'ND':1, 'IEA':1}


class Gonet(object):

    def __init__(self, geneOnly= True, filt=['other','unknown']\
    ,path='/Users/xiaoyiou/Documents/research/GO/ATH_GO_GOSLIM.txt'):
        self.raw = go_parser.goReader(path)
        for filter in filt:
            self.raw = self.raw[~self.raw['slim'].str.contains(filter)]
        if geneOnly:
            self.raw = self.raw[self.raw.kind == 'gene']
        self.__createNets()
        print "created go term lists"
        self.g2go, self.gos = self.__loadData()
        print "done with loading data"
        #self.dists = self.__calcDists()
        #self.go2id, self.id2go = None, None

    def __createNets(self, col='goid'):
        self.C = set(self.raw[self.raw['aspect'] == 'C'][col].unique().tolist())
        self.F = set(self.raw[self.raw['aspect'] == 'F'][col].unique().tolist())
        self.P = set(self.raw[self.raw['aspect'] == 'P'][col].unique().tolist())

    def __loadData(self, useslim=True):
        g2go = {}
        gos = {}
        for _, v in self.raw.iterrows():
            code, aspect, gene, desc, goid, kind, rel, slim = \
            v.code, v.aspect, v.gene, v.godesc, v.goid, v.kind, v.rel, v.slim
            if goid not in gos:
                gos[goid] = (aspect, rel, kind, slim, desc)
            if gene not in g2go:
                g2go[gene] = {}
            if goid not in g2go[gene]:
                g2go[gene][goid] = [code]
            else:
                g2go[gene][goid].append(code)
        return g2go, gos

    def __calcDists(self):
        dists = {}
        for lst in self.g2go.values():
            for i in xrange(len(lst)):
                for j in xrange(i+1, len(lst)):
                    key = frozenset({lst[i], lst[j]})
                    if key not in dists:
                        dists[key] = 1
                    else:
                        dists[key] += 1
        return dists

    def getDesc(self, glst):
        report = dict()
        for gene in glst:
            if gene in self.g2go:
                for goid in self.g2go[gene]:
                    if goid not in report:
                        report[goid] = 0
                    codes = [evidence[code] \
                             if (code in evidence) else 0 \
                             for code in self.g2go[gene][goid]]

                    report[goid] += sum(codes)
        res = []
        for x in sorted(report.items(), key = lambda x: x[1]):
            res.append((self.gos[x[0]], x[1]/float(len(glst))))
        return res[::-1]


    def drawNet(self, path, alpha=0):
        thresh = alpha * np.std(self.dists.values()) + np.mean(self.dists.values())
        g = gt.Graph(directed=False)
        go2id = {}
        id2go = {}
        id = 1
        for i, v in self.dists.items():
            if v >= thresh:
                a, b = i
                if a not in go2id:
                    go2id[a] = id
                    id2go[id] = a
                    id += 1
                if b not in go2id:
                    go2id[b] = id
                    id2go[id] = b
                    id += 1
                g.add_edge(go2id[a], go2id[b])

        pos = gt.sfdp_layout(g, multilevel = True)
        gt.graph_draw(g, pos, output_size= (1000,1000), output=path)
        self.go2id, self.id2go = go2id, id2go


    def report(self, glst, k=5):
        """

        Parameters
        ----------
        glst    : the gene list

        Returns : The top k items in each category
        -------

        """
        lst  = self.getDesc(glst)
        res = {'C':[], 'F':[], 'P':[]}
        N = k * 3
        i = 0
        while i < N:
            if i >= len(lst):
                break
            t, score = lst[i]
            a, kind, _, _, desc = t
            i += 1
            if len(res[a]) >= k:
                N += 1
                continue
            res[a].append((kind, desc, score))
        return res


    def printReport(self, glst, k=5, ):
        report = self.report(glst, k=k)
        for key in report:
            lst = report[key]
            if key == 'C':
                print '----------Components--------------'
            elif key == 'F':
                print '----------Function----------------'
            else:
                print '----------Process-----------------'

            for entry in lst:
                rel, desc, score = entry
                print "%s %s with score=%.3f"%(rel, desc, score)




