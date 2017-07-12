import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from itertools import cycle

linecycler= cycle(["-","--","-.",":","-v"])
methods = [
    'catworks',
    'svm',
    'GBN',
    'LogReg',
    'BAM'
    ]
mdata = np.array([np.mean(x, axis=0) for x in pmRatios])
sdata = np.array([np.std(x, axis=0) for x in pmRatios])

plt.figure()
xs = np.arange(0.1, 1.1, 0.1).tolist()
for i in range(mdata.shape[1]):
    mplot = plt.plot(xs, mdata[:,i], next(linecycler),label=methods[i], linewidth=2)
    dlower = np.clip(mdata[:,i] + sdata[:,i], 0, 1)
    dupper = np.clip(mdata[:,i] - sdata[:,i], 0, 1)
    plt.fill_between(xs,dlower,dupper,alpha=0.2, color=mplot[0].get_color())




plt.ylim([0.0, 1.05])
plt.xlabel('Pattern Ratio', fontsize=16)
plt.ylabel('Perfect Matching Ratios', fontsize=16)
plt.legend(fontsize=16)
plt.show()


