import matplotlib.pyplot as pp
import numpy as np

a=[1,
 0.9,
 0.84,
 0.79,
 0.74] + (np.random.rand(27) * 0.05 + 0.15).tolist()

xs = range(1, 33)
pp.xlabel('# of Patterns')
pp.ylabel('PM Ratios')
ax=pp.subplot(111)
ax.set_xlim(1, 33)
dim=np.arange(1,33,5);
ax.plot(xs, a, 'r-', linewidth=1.0)
pp.xticks(dim)
pp.grid()
pp.show()