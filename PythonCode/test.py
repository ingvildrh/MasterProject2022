'''
from datetime import datetime
from math import inf
import matplotlib.pyplot as plt
import numpy as np
import scipy 
from scipy import optimize
from scipy.sparse import csr_matrix
import scipy.sparse as sp
from scipy.optimize import minimize
'''

import datetime
import cProfile

from re import U
import time
from timeit import default_timer as timer
import numpy as np
arr = np.ones((100,1))
arr2 = np.zeros((100,1))
arr4 = np.zeros((100,1))
start = timer()
for i in range(100):
    arr2[i] = arr[i][0]*2
    arr4[i] = arr[i][0]*4
end= timer()
timetot = end-start
print(timetot)


'''
arr2 = arr*2
arr4 = arr*4


'''