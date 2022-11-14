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
import time
import datetime
import cProfile


def slow_add(a, b):
    time.sleep(0.5)
    return a+b

def fast_add(a, b):
    return a+b

prof = cProfile.Profile()

def main_func():
    arr = []
    prof.enable()
    for i in range(10):

        if i%2==0:
            arr.append(slow_add(i,i))
        else:
            arr.append(fast_add(i,i))
    prof.disable()
    return arr

main_func()