import numpy as np
import scipy
from scipy import optimize

def deb(A):
    print(A)
    print("dimention:")
    print(A.shape)

def termset(A,B,C,K,ymax,ymin,umax,umin):

    nx = np.shape(A)[0]
    nu = np.shape(B)[1]
    ny = np.shape(C)[0]

    big = 1000

    H0 = np.concatenate((np.identity(nx), -1*np.identity(nx)))
    h0 = big*np.ones((2*nx, 1))

    Hi  = np.concatenate((C, -1*C, -1*K, K))
    hi = np.concatenate((ymax, -1*ymin, np.array([[umax], [-1*umin]])))
    nci = hi.shape[0]

    ni = 10
    nstep = -1

    xlist = []
    fvallist = []
    while(ni>0):
        ni = 0
        nstep = nstep+1
        if (nstep==0):
            Dyn = np.identity(nx)
        else:
            Dyn = (np.subtract(A, (B.dot(K)))).dot(Dyn)
        for ik in range(nci):
            f = (np.matrix((Hi[ik, :]).dot(Dyn)))
            f = np.transpose(f)
            OPT = optimize.linprog(-1*f, A_ub = H0, b_ub=h0) 
            x = OPT.x #the optimal x, dont think this is used any further, might remove it
            fval = OPT.fun #objective func value at minimum
            xlist.append(x)
            fvallist.append(fval)
            if (-1*fval > hi[ik]):
                ni = ni +1
                H0 = np.concatenate((H0, np.transpose(f)))
                print(hi[ik])
                h0 = np.append(h0, hi[ik])
                
            
    
    nc0 = np.shape(H0)[0] #generelt litt usikker på disse med shape

    H1 = H0[2*nx+1:nc0, :]
    h1 = h0[2*nx+1:nc0]

    for ik in range(2*nx):
        f = np.transpose(H0[ik, :])
        OPT = scipy.optimize.linprog(-1*f, H0, h0)
        fval = OPT.fun
        if(-1*fval > h0[ik]):
            H1 = np.concatenate(H1, np.transpose(f))
            h1 = np.concatenate(h1, h0[ik])
    
    Hn = H1
    hn = h1
    return Hn, hn







