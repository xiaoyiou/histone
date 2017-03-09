# This script is a test for hetero association learning using non-zero labels

from eonetwork import mapper,eval

learner = mapper.Assoc(hopfield=True, keepz=False)
A, B = learner.loaddata(all, go, 'F')

# removing all-zero entries
inds = [x[0] for x in(B==1).any(axis=1).tolist()]
A = A[inds, :]
B = B[inds, :]

learner.train(A, B)