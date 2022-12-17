from re import U
import numpy as np
from MPC_setup import * 
import math as math
import time
import daug
import datetime
import scipy.io
import cProfile

from timeit import default_timer as timer

#Time steps for the MPC
tend = 100

#code added to use line profiler in the lineproftest.py script

#def main():
if (MODEL == 1):
    xinit = np.matrix([[5], [-2]])
    tosave1 = np.zeros((1, tend))
if (MODEL == 2):
    xinit = np.matrix([[25.5724],[25.3546],[9.7892],[0.2448]])
    tosave2 = np.zeros((1, tend))

x0 = xinit

#Allocating memory for filling in as the computations finish
xsave = np.zeros((nx, tend+1))
usave = np.zeros((nu, tend))

feasflag = 1

xsave[:,0:1]= xinit

nc = np.shape(Gz)[0]
Qreset = np.identity(nc) #improvement 3: making the Qreset matrix once, outside of the loop
Hinv = np.linalg.inv(Hess) 
GzHinv = Gz.dot(Hinv)
GzHinvGz = GzHinv.dot(np.transpose(Gz))
IGIs = np.subtract(Qreset, GzHinvGz)
hif = Hinv@f0

HGz = -Hinv@np.transpose(Gz) 
HGzu = HGz[0:nu, :]
hifu = -1*hif[0:nu, :]

while_iterations1 = [] 
while_iterations2 = []

tot = timer()

for i in range(tend):
    start = timer()
    feasflag = False
    solved = False
    ix = 0

    #Qmat0i = np.identity(nc) before improvement 3
    Qmat0i = Qreset #improvement 3: using the matrix made outside the loop
    actset = np.zeros((nc,1))

    y0 = np.subtract(-Sz.dot(x0), Wz)
    while (not (solved)):
        ix = ix +1
        if (ix == 1): 
            y = y0
        else:
            y = np.subtract((y0), (y0[iz].item())*vAd)
            y0 = y
        lam = np.multiply((y),actset)
        #i1 = min(lam) before improvement 1
        i1=lam.min() #improvement 1: using a NumPy function with vectorization
        i1z = np.argmin(lam)
        if (i1>=0):
            i1 = None
        if (i1 != None):
            iz = i1z
            actset[iz] = 0 
            qc = 1
        else:
            #i2 = max(y-lam) before improvement 1
            y_lam = y-lam
            i2 = (y_lam).max() #improvement 1: using a NumPy function with vectorization
            i2z = np.argmax(y_lam)
            if (i2 <= 0):
                i2= None
            if (i2 != None): 
                iz = i2z
                actset[iz] = 1
                qc = -1
            else:
                iz = None
        
        if (iz != None):
            qu = IGIs[:, iz]
            vA = np.transpose(np.matrix(Qmat0i@(qu)))
            qdiv = qc+vA[iz] 
            
            if ((abs(qdiv) < 1*math.e**(-13)) or (abs(qdiv) > 1*math.e**(12))):

                feasflag = 0
                print("Infeasible problem detected")

                feasflag = False

                break

            #vAd = np.multiply((1/qdiv),(vA)) before improvement 1
            vAd = np.divide((vA), qdiv) #improvement 1
            
            vAdQmat0i = vAd@np.matrix(Qmat0i[iz,:])
            Qmat1i = np.subtract(Qmat0i,  vAdQmat0i)

            Qmat0i = Qmat1i

        else:
            solved = True
            feasflag = True
    
    
    end = timer()
    tk = end-start
    if (MODEL==1):
        while_iterations1.append(ix)
    if (MODEL==2):
        while_iterations2.append(ix)
    if (feasflag == 0):
        break
    
    u = HGzu@lam+hifu@x0
    x1 = A@x0+B@u

    x0 = x1
    xsave[:, i+1] = np.transpose(x0)
    usave[:, i] = np.transpose(u)

    if (MODEL ==1):
        tosave1[0, i] = tk 
    if (MODEL ==2):
        tosave2[0, i] = tk 

end_tot = timer()
total_time = end_tot-tot

#code added to use line profiler in the lineproftest.py script
#main()


