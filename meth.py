# This script is used to visualize the DNA methylation
import pandas as pd
import numpy as np

meth = pd.read_csv('dnaMethy.csv', header = 0, index_col=9)

meth_labeled = pd.merge(meth,young.labels, right_index=True, left_index=True)

cut_off = .5

meth_labeled = meth_labeled.replace(['nd', 'Inf', \
                          '-Inf' ], np.NaN)


meth_labeled = meth_labeled.convert_objects(convert_numeric=True).dropna()



print meth_labeled[(meth_labeled['pvalue'] <= .05)&(meth_labeled['flowerN'] == 1)\
             &(meth_labeled['defense'] == 1)]