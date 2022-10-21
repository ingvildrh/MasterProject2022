import numpy as np
from MPC_run import *
import matplotlib.pyplot as plt

exe_time = np.matrix(tosave)
mlist = []
#as matlab returned NaN when reading the first and last number, I add a buffer in both ends
mlist.append("buff1")

for i in range(100):
   mlist.append(exe_time[0, i])

mlist.append("buff2")

if (MODEL == 1):
    file = open("execution_timeM1.txt", "w")
    file.write(str((mlist)))

if (MODEL == 2):
    file = open("execution_timeM2.txt", "w")
    file.write(str((mlist)))