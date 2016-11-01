# This file plots the correlation analysis of pairs of

# plotting pScore vs backgroun ratios

import numpy as np
import matplotlib.pyplot as plt


px = young.ratios['all']

fig =plt.figure()
fig.suptitle("Score vs Global ratios", fontsize=16)


ax= plt.subplot("611")
ax.set_title("Stress")
py = young.pScore['stress']
ax.scatter(px,py)
a,b = np.polyfit(px,py,1)
plt.plot(px,a*px+b,'-')

ax= plt.subplot("612")
ax.set_title("Development")
py = young.pScore['develop']
ax.scatter(px,py)
a,b = np.polyfit(px,py,1)
plt.plot(px,a*px+b,'-')
ax= plt.subplot("621")

ax.set_title("Stimulus")
py = young.pScore['stimulus']
ax.scatter(px,py)
a,b = np.polyfit(px,py,1)
plt.plot(px,a*px+b,'-')

ax= plt.subplot("622")
ax.set_title("Flowring")
py = young.pScore['flowerN']
ax.scatter(px,py)
a,b = np.polyfit(px,py,1)
plt.plot(px,a*px+b,'-')

ax= plt.subplot("631")
ax.set_title("defense")
py = young.pScore['defense']
ax.scatter(px,py)
a,b = np.polyfit(px,py,1)
plt.plot(px,a*px+b,'-')




plt.show()

