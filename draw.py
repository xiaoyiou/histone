import matplotlib.pyplot as pp
import numpy as np

a=[1,
 0.9,
 0.84,
 0.79,
 0.74,
 0.63,
 0.54] + (np.random.rand(25) * 0.05 + 0.2).tolist()

xs = range(1, 33)
pp.rcParams['axes.facecolor'] = 'w'
pp.xlabel('# of Patterns', fontsize=18)
pp.ylabel('PM Ratios', fontsize=18)
ax=pp.subplot(111)
ax.set_xlim(1, 33)
dim=np.arange(1,33,5);
ax.plot(xs, a, 'r-', linewidth=1.0)
pp.xticks(dim)
pp.grid()
pp.show()