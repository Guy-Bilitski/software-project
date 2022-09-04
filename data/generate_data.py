from sklearn.datasets import make_blobs
import numpy as np

def sym(n):
    a = 10*np.random.rand(n, n)
    m = np.tril(a) + np.tril(a, -1).T
    np.savetxt("data/mytxt.txt", m, fmt='%1.3f', delimiter=",")

def data(samples=10, features=5, centers=3):
    #a = np.random.rand(n, m)
    X, y = make_blobs(n_samples=samples, centers=centers, n_features=features, random_state=0)
    np.savetxt("data/mytxt.txt", X, fmt='%1.3f', delimiter=",")
