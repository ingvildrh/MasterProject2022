import numpy as np
from MPC_run import *
import matplotlib.pyplot as plt

if (MODEL ==1):
    exe_time = np.matrix(tosave1)
    iterations = while_iterations1
if (MODEL ==2):
    exe_time = np.matrix(tosave2) 
    iterations = while_iterations2

time_list = []
iteration_list = []
#as matlab returned NaN when reading the first and last number, I add a buffer in both ends
time_list.append("buff1")
iteration_list.append("buff1")

for i in range(100):
   time_list.append(exe_time[0, i])
   iteration_list.append(iterations[i])

time_list.append("buff2")
iteration_list.append("buff1")

if (MODEL == 1):
    file = open("..//Matlab_kode//execution_timeM1.txt", "w")
    file.write(str((time_list)))

    file = open("..//Matlab_kode//while_iterations1.txt", "w")
    file.write(str((iteration_list)))

    print(iteration_list)
    print(time_list)

if (MODEL == 2):
    file = open("..//Matlab_kode//execution_timeM2.txt", "w")
    file.write(str((time_list)))

    file = open("..//Matlab_kode//while_iterations2.txt", "w")
    file.write(str((iteration_list)))

    print(iteration_list)
    print(time_list)