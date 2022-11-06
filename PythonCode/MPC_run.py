from re import U
import numpy as np
from MPC_setup import * 
import math as math
import time
import daug
import datetime
import scipy.io

tend = 100

if (MODEL == 1):
    xinit = np.matrix([[5], [-2]])
    tosave1 = np.zeros((1, tend))
if (MODEL == 2):
    xinit = np.matrix([[25.5724],[25.3546],[9.7892],[0.2448]])
    tosave2 = np.zeros((1, tend))

x0 = xinit

xsave = np.zeros((nx, tend+1))
usave = np.zeros((nu, tend))
zsave = np.zeros((nu, tend))

feasflag = 1

xsave[:,0:1]= xinit

nc = np.shape(Gz)[0]
Hinv = np.linalg.inv(Hess) #stemmer
IGi = np.subtract(np.identity(nc), (Gz.dot(Hinv)).dot(np.transpose(Gz)))
IGIs = (IGi) #stemmer, but not sparse
hif = Hinv@f0 #stemmer

HGz = -Hinv@np.transpose(Gz) #stemmer
HGzu = HGz[0:nu, :]
hifu = -hif[0:nu, :]

actset0 = np.zeros((nc,1))
Qmat0i = np.identity(nc) #denne er overflødig siden vi setter dne i for loopen

while_iterations1 = [] #for counting the while iterations in each step
while_iterations2 = []

for i in range(tend):
    start = time.time_ns()    
    
    feasflag = False

    actset = actset0
    solved = False
    ix = 0

    Qmat0i = np.identity(nc)
    
    y0 = np.subtract(-Sz.dot(x0), Wz)

    
    while (not (solved)):
        
        ix = ix +1
        
        if (ix == 1): 
            y = y0
        else:
            y = np.subtract((y0), (y0[iz].item())*vAd)
            y0 = y
        
        lam = np.multiply((y),actset) #elementvis?
        i1 = min(lam)
        i1z = np.argmin(lam)

        if (i1>=0):
            i1 = []

        if (i1):
            iz = i1z
            actset[iz] = 0 
            qc = 1
        else:
            i2 = max(y-lam)
            i2z = np.argmax(y-lam)
            if (i2 <= 0):
                i2= []
            if (i2): 
                iz = i2z
                actset[iz] = 1
                qc = -1
            else:
                iz = []
        
        if (iz):
            qu = IGIs[:, iz]
            vA = np.transpose(np.matrix(Qmat0i@(qu)))
            qdiv = qc+vA[iz] 
            
            if ((abs(qdiv) < 1*math.e**(-13)) or (abs(qdiv) > 1*math.e**(12))):

                feasflag = 0
                print("Infeasible problem detected")

                feasflag = False

                break

            vAd = np.multiply((1/qdiv),(vA))
            #Qmat1i = np.subtract(Qmat0i, vAd@np.matrix(Qmat0i[iz,:]))
            o=vAd@np.matrix(Qmat0i[iz,:])
            Qmat1i = np.subtract(Qmat0i, vAd@np.matrix(Qmat0i[iz,:]))
            #Qmat1i = Qmat0i-vAd@np.matrix(Qmat0i[iz,:])
            Qmat0i = Qmat1i

        else:
            solved = True
            feasflag = True
    
    
    end = time.time_ns()
    tk = end-start
    if (MODEL==1):
        while_iterations1.append(ix)
    if (MODEL==2):
        while_iterations2.append(ix)
    if (feasflag == 0):
        break
    
    #lam = np.multiply((actset),(y)) #element wise?

    u = HGzu@lam+hifu@x0
    x1 = A@x0+B@u

    x0 = x1
    xsave[:, i+1] = np.transpose(x0)
    usave[:, i] = np.transpose(u)

    if (MODEL ==1):
        tosave1[0, i] = tk #this is in nano seconds

    if (MODEL ==2):
        tosave2[0, i] = tk #this is in nano seconds
        

print("done")




