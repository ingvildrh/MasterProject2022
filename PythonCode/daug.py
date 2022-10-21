import numpy as np

def daug(A, B):
    Am = np.shape(A)[0]
    An = np.shape(A)[1]
    Bm = np.shape(B)[0]
    Bn = np.shape(B)[1]
    C = np.zeros((Am+Bm, An+Bn))
    C[0:Am, 0:An]=A
    C[Am:Am+Bm, An:An+Bn]=B
    return C