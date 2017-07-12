"""
This module is used to test different methods for multi-label classification problems
"""


# Baseline, transform this into a single class problem by creating
# artificial labels

from sklearn.naive_bayes import GaussianNB as GNB
from skmultilearn.problem_transform import BinaryRelevance
from sklearn.linear_model import LogisticRegression as LR
from sklearn.svm import SVC
from sklearn.metrics import hamming_loss
from sklearn.ensemble import BaggingClassifier, RandomForestClassifier
import numpy as np
import catworks.simulate as sim
import catworks.catworks as cts
import catworks.bam as bam
import timeit
import gc
def getpm(y, py):
    return np.sum(np.all(np.equal(y, py), axis=1)) / 1.0 / len(py)




l, r, n, pr, pn, rn = 15, 15, 20000, 1, 150, 0.15
prs = np.arange(1.0, 1.01, 0.1)
names = ['catworks', 'SVM', 'Naive Bayeisan',
         'Logistic Regression', 'BAM']

pmRatios = []
hmLoss = []
rtimes = []


for pr in np.arange(0.1, 1.01, 0.1):
    times, pms, hms = [], [], []
    for _ in xrange(10):
        X, Y, px, py = sim.createData(l, r, n, 0, pr, [1.0/pn]*pn, rn)
        py = np.array(py)
        time = []
        pm = []
        hm = []

        # catworks
        start = timeit.default_timer()
        cat = cts.Acat(l, r)
        cat.feed(X, Y)
        y = cat.getRes(cat.meow(px, thresh=2, k=3))
        end = timeit.default_timer()
        time.append(end - start)
        hm.append(hamming_loss(y, py))
        pm.append(getpm(y, py))

        # SVM
        n_estimators = 5
        clf = BaggingClassifier(SVC(kernel='linear', probability=True), max_samples=1.0 / n_estimators,
                      n_estimators=n_estimators)
        start = timeit.default_timer()
        model = BinaryRelevance(clf, require_dense =True)
        model.fit(X, Y)
        y = model.predict(px).toarray()
        end = timeit.default_timer()
        time.append(end - start)
        hm.append(hamming_loss(y, py))
        pm.append(getpm(y, py))

        # GBN
        start = timeit.default_timer()
        model = BinaryRelevance(GNB(), require_dense=True)
        model.fit(X, Y)
        y = model.predict(px).toarray()
        end = timeit.default_timer()
        time.append(end - start)
        hm.append(hamming_loss(y, py))
        pm.append(getpm(y, py))

        # Logistirc Regression
        start = timeit.default_timer()
        model = BinaryRelevance(LR(), require_dense=True)
        model.fit(X, Y)
        y = model.predict(px).toarray()
        end = timeit.default_timer()
        time.append(end - start)
        hm.append(hamming_loss(y, py))
        pm.append(getpm(y, py))
     #    BAM
        start = timeit.default_timer()
        bm = bam.BAM(bam.trans_data(X, Y))
        y = bm.predict(px)
        end = timeit.default_timer()
        time.append(end - start)
        hm.append(hamming_loss(y, py))
        pm.append(getpm(y, py))
        times.append(time)
        pms.append(pm)
        hms.append(hm)
        gc.collect()

    pmRatios.append(np.array(pms))
    hmLoss.append(np.array(hms))
    rtimes.append(np.array(time))
