"""
hofield module is for hetero/auto association learning for patterns
"""
from pylab import imshow, cm, show
import numpy as np




class Hopfield(object):
    def __init__(self, keepZero=False):
        self.r, self.c = 0, 0
        self.keepz = keepZero
        self.W = None

    def train(self, data):
        self.r, self.c = data.shape
        self.W = np.zeros((self.c, self.c))
        for p in data:
            self.W = self.W + np.outer(p, p)
        self.W[np.diag_indices(self.c)] = 0
        self.W = self.W/self.r

    def recall(self, patterns, steps= 20):
        if self.keepz:
            sgn = np.vectorize(lambda x: -1 if x < 0
                                else (0 if x == 0 else +1))
        else:
            sgn = np.vectorize(lambda x: -1 if x < 0 else +1)
        for _ in xrange(steps):
            cur = sgn(np.dot(patterns, self.W))
            if (cur == patterns).all():
                break
            patterns = cur
        return patterns


class HeteroAssoc(object):

    def __init__(self, keepZero = False):
        self.i, self.o = 0, 0
        self.keepz = keepZero
        self.W = None

    def train(self, I, D):
        """
        Parameters
        ----------
        I : The matrix of n-by-N where n is dimension of input
        D : The matrix of m-by-N where m is the dimension of the output

        Returns: Nothing. Just update the weight matrix
        -------
        """
        if I.shape[1] != D.shape[1]:
            print "unmatching dimensions",I.shape,D.shape
            return
        i, d = np.mat(I), np.mat(D)
        self.W = d * i.T * np.linalg.pinv(i * i.T)

    def recall(self, I):
        """

        Parameters
        ----------
        I : the n-by-N where n is the dimesion of inputs

        Returns
        -------

        """
        return  self.W * np.mat(I)

if __name__ == '__main__':
    A = """
    .XXX.
    X...X
    XXXXX
    X...X
    X...X
    """

    Z = """
    XXXXX
    ...X.
    ..X..
    .X...
    XXXXX
    """
    def to_pattern(letter):
        return np.array([+1 if c == 'X' else -1 for c in
                         letter.replace(' ','').replace('\n', '')])


    def display(pattern):
        imshow(pattern.reshape((5, 5)), cmap=cm.binary, interpolation='nearest')
        show()

    aa, zz= to_pattern(A), to_pattern(Z)

    display(aa);
    display(zz)