import numpy as np
import matplotlib.pyplot as plt

SLOW = 0 

if (SLOW):
    from MPC_run_slow import *
else:
    from MPC_run import *

PYTHONPLOTS = 1




if (MODEL ==1):
    exe_time = np.matrix(tosave1)
    iterations = while_iterations1
if (MODEL ==2):
    exe_time = np.matrix(tosave2) 
    iterations = while_iterations2

average_iterations = sum(iterations)/len(iterations)
print("average iterations: ")
print(average_iterations)


#Plotting runtime for Python algorithm 
if (PYTHONPLOTS):
    if (MODEL == 1):
        plt.figure(3)
        plt.plot(np.arange(tend), tosave1[0,:], color="green", label="Time")
        plt.title("Runtime Python Model " + str(MODEL))
        plt.xlabel("time step")
        plt.ylabel("Runtime")
        plt.grid()

    if (MODEL == 2):
        plt.figure(3)
        plt.plot(np.arange(tend), tosave2[0,:], color="green", label="Time")
        plt.title("Runtime Python Model " + str(MODEL))
        plt.xlabel("time step")
        plt.ylabel("Runtime")
        plt.grid()

plt.figure(4)
plt.plot(np.arange(tend), iterations, color="blue", label="Time")
plt.title("While iterations Python  " + str(MODEL))
plt.grid()
plt.xlabel("timestep")
plt.ylabel("number of iterations") 


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


#Just a plot similar to the one in matlab so that we can compare states and inputs with matlab
if (PYTHONPLOTS):
    if (MODEL == 1):
        states = plt.figure(1)
        plt.plot(np.arange(tend+1),xsave[0], color="red", label="x1")
        plt.plot(np.arange(tend+1), xsave[1], color="blue", label="x2")
        plt.title("State trajectories")
        plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        plt.xlim([0,100])
        plt.ylim([-3,5])
        plt.yticks([-3,-2,-1,0,1,2,3,4,5])
        plt.legend(loc="lower center",fontsize=10,ncol=2)
        plt.xlabel("time step")
        plt.ylabel("x")
        plt.grid()
        

        input = plt.figure(2)
        plt.plot(np.arange(tend), usave[0], color="blue", label="u")
        plt.title("System input")
        plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        plt.xlim([-5, 100])
        plt.ylim([-0.6,1.1])
        plt.yticks([-0.6,-0.4,-0.2,-0,0.2,0.4,0.6,0.8,1])
        plt.xlabel("time step")
        plt.ylabel("u")
        plt.legend("u")
        plt.grid()
        plt.show()

    if (MODEL == 2):
        ym = C@xsave 
        output_trajectories = plt.figure(1)
        plt.plot(np.arange(tend+1), ym[0], color="blue", label="ym1")
        plt.plot(np.arange(tend+1), ym[1], color="red", label="ym2")
        plt.title("Output trajectories")
        plt.xlim([0,100])
        plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        plt.ylim([-1.5, 1])
        plt.yticks([-1.5,-1,-0.5,0,0.5,1])
        plt.xlabel("time step")
        plt.ylabel("ym")
        plt.legend(loc="lower center",fontsize=10,ncol=2)
        plt.grid()
        
        ufig = plt.figure(2)
        plt.plot(np.arange(tend), usave[0], color="blue", label="u1")
        plt.plot(np.arange(tend), usave[1], color="red", label="u2")
        plt.title("Input trajectories")
        plt.xlim([0,100])
        plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        plt.ylim([-0.7,0])
        plt.xlabel("time step")
        plt.legend(loc="lower center",fontsize=10,ncol=2)
        plt.ylabel("u")
        plt.grid()
        plt.show()

