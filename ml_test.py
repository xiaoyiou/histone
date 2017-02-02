"""
This module is used to test different methods for multi-label classification problems
"""


# Baseline, transform this into a single class problem by creating
# artificial labels

from sklearn.naive_bayes import GaussianNB
from skmultilearn.problem_transform import BinaryRelevance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import hamming_loss
from sklearn.model_selection import StratifiedKFold as KFold
from random import shuffle
# Create dataset with reasonable size of unlabelled data
Y = young.labels[young.labels.columns[1:]]
#Y = Y[(Y.T !=0).any()]
Y = Y.ix[Y.sum(axis = 1) > 1]
XX = X.ix[Y.index]

precs = []
indexs = range(XX.shape[0])

for _ in xrange(5):
    index_shuffle = shuffle(indexs)
    XX = XX.iloc[indexs,:]
    Y = Y.iloc[indexs,:]
    train = XX.shape[0] / 4
    model = BinaryRelevance(GaussianNB(), require_dense = True)
    model.fit(XX.iloc[0:train,:],Y.iloc[0:train,:])
    prediction = model.predict(XX.iloc[train:,:])
    precs.append(hamming_loss(prediction.toarray(),Y.iloc[train:,:]))


