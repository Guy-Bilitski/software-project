import numpy as np
from math import sqrt
import mykmeanssp

def max_elem(A: np.array):
    d = {}
    n = A.shape[0]
    maxi = float('-inf')
    ii = -1
    jj = -1
    for i in range(n):
        for j in range(i+1,n):
            if abs(A[i,j]) > maxi:
                ii = i
                jj = j
                maxi = abs(A[i,j])
    return ii, jj, maxi


def s_and_c(i,j,maxi, A):
    theta = (A[j,j]-A[i,i])/(2*A[i,j])
    sign = 1 if theta >=0 else -1
    t = sign / (abs(theta) + sqrt(theta**2 + 1))
    c = 1/sqrt(t**2 + 1)
    s = t*c
    return s,c
            
def off(A):
    n = A.shape[0]
    of = 0
    for i in range(n):
        for j in range(n):
            if i!=j:
                of += A[i,j]**2
    return of


def create_p(s,c,i,j,maxi, n):
    P = np.eye(n)
    P[i,j] = s
    P[j,i] = -s
    P[i,i] = c
    P[j,j] = c
    return P



def jacobi(A):
    n = A.shape[0]
    V = np.eye(n)
    
    for i in range(100):
        old_off = off(A)
        i,j,maxi = max_elem(A)
        s,c = s_and_c(i,j,maxi,A)
        P = create_p(s,c,i,j,maxi,n)
        V = V@P
        A = (P.T) @ A @ P
        new_off = off(A)
        if old_off - new_off < 0.00001:
            break
    return A,V

    
    
def eigen_gap_h(A,V):
    n = A.shape[0]
    evals = sorted([(A[i,i], i) for i in range(n)])
    evals = evals[::-1]
    maxi = -1
    k = -1
    for i in range(n-1):
        if evals[i][0] - evals[i+1][0] > maxi:
            maxi = evals[i][0] - evals[i+1][0]
            k = i
    return k+1
            



def load_initial_A_for_jacobi(path='testfiles/spk_2.txt'):
    data_points = np.genfromtxt(fname=path, dtype=float, delimiter=',')
    A = mykmeanssp.lnorm(data_points.tolist())
    return np.array(A)

A = load_initial_A_for_jacobi()
A,V = jacobi(A)
print(A.shape)
k=eigen_gap_h(A,V)
print(f'PYTHON:\tk={k}')



#np.savetxt("data/Tpy.txt", V, fmt='%1.3f', delimiter=",")