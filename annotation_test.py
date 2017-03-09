from __future__ import print_function
import gonet
import matplotlib.pyplot as plt
sorted_all = sorted(all.saxGrps.items(), key=lambda x: len(x[1]))[::-1]
go = gonet.Gonet(geneOnly=False)
min_size = 10
# Printing the patterns
path = 'Annotation_code.txt'
xx = []
yy1 = []
yy2 = []
yy3 = []
with open(path,'wb') as f:
    for pattern, glst in sorted_all:
        if len(glst) < min_size:
            break
        xx.append(len(glst))
        print("##############################################",file = f)
        print([int(i) for i in pattern], len(glst), file=f)
        print("##############################################", file=f)
        report = go.report(glst, k=6)
        for key in report:
            lst = report[key]
            if key == 'C':
                print ('----------Components--------------',file=f)
            elif key == 'F':
                print ('----------Function----------------', file=f)
            else:
                print ('----------Process-----------------', file=f)

            for entry in lst:
                rel, desc, score = entry
                print ("%s %s with score=%.3f" % (rel, desc, score), file=f)
        yy1.append(report['C'][0][-1])
        yy2.append(report['F'][0][-1])
        yy3.append(report['P'][0][-1])
        print('\n' ,file=f)
    f.close()

f, axarr = plt.subplots(2,5)
for i in xrange(10):
    N = len(all.saxnets[i].cluster_centers)
    for j in xrange(N):
        s = all.saxnets[i].cluster_centers[j]
        axarr[i/5, i%5].plot(s, label=j)
        axarr[i/5, i%5].legend()
    axarr[i/5,i%5].set_title(mod_names[i])