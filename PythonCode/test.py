from datetime import datetime
from math import inf
import matplotlib.pyplot as plt
import numpy as np
import scipy 
from scipy import optimize
from scipy.sparse import csr_matrix
import scipy.sparse as sp
from scipy.optimize import minimize
import time
import datetime

'''
ANS = scipy.optimize.linprog(
    c = [-1, -2], 
    A_ub=[[1, 1]], 
    b_ub=[6],
    bounds=(1, 5),
    method='highs'
)
'''
start = time.process_time()
start1 = datetime.datetime.now()

i = 1000000
while (i != 0):
    i = i-1


#tk = time.process_time() - start 
tk1 = datetime.datetime.now()-start1

#print(tk)
print(tk1)



# f= np.array([[0, -1]])
# H0 = np.array([[1,0], [0,1 ], [-1,0], [0,-1], [1,0], [0,1]])
# h0=np.array([[1000], [1000], [1000],[1000],[5],[5]])
# h01=np.array([[1000, 1000, 1000,1000,5,5]])
# OPT = optimize.linprog(-1*f, A_ub = H0, b_ub=h0) 
# x = OPT.x #the optimal x, dont think this is used any further, might remove it
# fval = OPT.fun 
# print(x)
# print(fval)

# fun = lambda x: -1*x[1]
# cons = ({'type': 'ineq', 'fun': lambda x: 1000-x[0]},
# {'type': 'ineq', 'fun': lambda x: 1000-x[1]}, 
# {'type': 'ineq', 'fun': lambda x: 1000+x[0]}, 
# {'type': 'ineq', 'fun': lambda x: 1000+x[1]}, 
# {'type': 'ineq', 'fun': lambda x: 5-x[0]},
# {'type': 'ineq', 'fun': lambda x: 5-x[1]})
# sol = minimize(fun, (0,0), method = 'SLSQP', constraints=cons)
# print(sol.x)








# Q = np.array([[-1], [0]])
# A = np.array([[1,0], [0,1 ], [-1,0], [0,-1], [1,0], [0,1]])
# b = np.array([[1000], [1000], [1000],[1000],[5],[5]])
# m = gp.Model()
# x = m.addMVar(2, ub=1.0)

# m.setMObjective(Q, b, xQ_L=None, xQ_R=None, xc=None, sense=None)
# m.addMCinstrs(A, x, sense, b)
# m.optimize()

# print(x.X)
# print('Obj: %g' % m.ObjVal)


# OPT = optimize.linprog(-1*f, A_ub = H0, b_ub=h0, bounds=(None, 0), method='highs-ds') 
# print(OPT.fun)
# print(OPT.x)
# print(OPT.success)

# def daug(A, B):
#     Am = np.shape(A)[0]
#     An = np.shape(A)[1]
#     Bm = np.shape(B)[0]
#     Bn = np.shape(B)[1]
#     C = np.zeros((Am+Bm, An+Bn))
#     C[0:Am, 0:An]=A
#     C[Am:Am+Bm, An:An+Bn]=B
#     return C


# a = np.array([[9, 2, 3], [1, 2, 3], [2,4,5]])
# b = np.matrix([[2, 2, 3], [1, 2, 3], [5,8,9]])
# d = np.array([[2, 2, 3], [2, 3, 4]])
# e = np.zeros((2,5))
# a[0:3, 0] = np.array(([1,4,5]))
# print(a)
# print("e")
# print(e)
# g = np.array([1, 1])
# print(g)
# e[:,1] = g
# print(e)
# A = sp.eye(3)
# F = e.dot(A)




