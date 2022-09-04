import numpy as np
from math import sqrt
from numpy import linalg as LA
import data.generate_data as gd

path = 'data/mytxt.txt'
gd.data()
X = np.genfromtxt(fname=path, dtype=float, delimiter=',')
print(X.shape)
gaussian_RBF = lambda x1, x2: np.exp(-LA.norm(x1-x2)/2)

def print_matrix(matrix):
        for row in matrix:
            line = []
            for value in row:
                line.append('%.4f,' % value)
            if len(line) > 0:
                line[-1] = line[-1][:-1]
                print("".join(line))
                
def py_wam(X):
    n = X.shape[0]
    W = [[None for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j:
                W[i][j] = 0
            else:
                W[i][j] = gaussian_RBF(X[i],X[j])
    print_matrix(W)
    return W

def py_ddg(X):
    W = py_wam(X)
    n = X.shape[0]
    D = [[0 for j in range(n)] for i in range(n)]
    for i, row in enumerate(W):
        D[i][i] = sum(row)
    return np.array(D)

def py_lnorm(X):
    D = py_ddg(X)
    W = py_wam(X)
    n = X.shape[0]
    for i in range(D.shape[0]):
        D[i,i] = 1 / (sqrt(D[i,i]))
    return np.eye(n) - D@W@D

    

print_matrix(py_lnorm(X))

# print()
# print(mykmeanssp.wam(X.tolist()))
# print()
# print(mykmeanssp.ddg(X.tolist()))
# print()
# print(mykmeanssp.lnorm(X.tolist()))