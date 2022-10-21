import numpy as np
from MPC_run import *
import matplotlib.pyplot as plt

if (MODEL ==1):
    exe_time = np.matrix(tosave1)
if (MODEL ==2):
    exe_time = np.matrix(tosave2) 


time_list = []

#as matlab returned NaN when reading the first and last number, I add a buffer in both ends
time_list.append("buff1")

for i in range(100):
   time_list.append(exe_time[0, i])

time_list.append("buff2")


if (MODEL == 1):
    file = open("..//Matlab_kode//execution_timeM1.txt", "w")
    file.write(str((time_list)))
    print(time_list)

if (MODEL == 2):
    file = open("..//Matlab_kode//execution_timeM2.txt", "w")
    file.write(str((time_list)))
    print(time_list)