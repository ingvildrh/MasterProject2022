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
    xinit = np.matrix([[-6], [1.1]])
if (MODEL == 2):
    xinit = np.matrix([[25.5724],[25.3546],[9.7892],[0.2448]])

x0 = xinit

xsave = np.zeros((nx, tend+1))
usave = np.zeros((nu, tend))
zsave = np.zeros((nu, tend))
tosave = np.zeros((1, tend))
feasflag = 1

xsave[:,0:1]= xinit

f = np.zeros((n*nu, 1))

nc = np.shape(Gz)[0]
Hinv = np.linalg.inv(Hess) #stemmer
IGi = np.subtract(np.identity(nc), (Gz.dot(Hinv)).dot(np.transpose(Gz)))
IGIs = (IGi) #stemmer, but not sparse
hif = Hinv@f0 #stemmer

HGz = -Hinv@np.transpose(Gz) #stemmer
HGzu = HGz[0:nu, :]
hifu = -hif[0:nu, :]

actsets = np.zeros((nc,1))
Qsp0 = np.identity(nc)



for i in range(tend):
    start = datetime.datetime.now()
    feasflag = False

    actset = actsets
    solved = False
    izold = 0
    ix = 0

    Qmat0i = Qsp0
    ylist = []
    vAdlist = []

    y0 = np.subtract(-Sz.dot(x0), Wz)
    while (not (solved)):

        ix = ix +1
        
        if (ix == 1): 
            y = y0
        else:
            q = (y0[iz].item())*vAd
            y01 = np.squeeze(y0)
            y = np.subtract((y0), (y0[iz].item())*vAd)
            y0 = y
        
        lam = np.multiply((y),actset) #elementvis?
        i1 = min(lam)
        i1z = np.argmin(lam)

        if (i1>=-10*np.finfo(float).eps):
            i1 = []

        if (i1):
            iz = i1z
            actset[iz] = 0 
            qc = 1
        else:
            O = y-lam 
            O2 = np.subtract(y, lam)
            i2 = max(y-lam)
            i2z = np.argmax(y-lam)
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
            vAd = np.multiply((1/qdiv),(vA))
            if ((abs(qdiv) < 1*math.e**(-13)) or (abs(qdiv) > 1*math.e**(12))):

                feasflag = 0
                print("Infeasible problem detected")

                feasflag = False

                break

            #Qmat1i = np.subtract(Qmat0i, vAd@np.matrix(Qmat0i[iz,:]))
            Qmat1i = np.subtract(Qmat0i, -vAd@np.matrix(Qmat0i[iz,:]))
            Qmat0i = Qmat1i

        else:
            solved = True
            feasflag = True
    tk = datetime.datetime.now() - start 
    if (feasflag == 0):
        break
    
    lam = np.multiply((actset),(y)) #element wise?

    u = HGzu@lam+hifu@x0
    x1 = A@x0+B@u

    x0 = x1
    xsave[:, i+1] = np.transpose(x0)
    usave[:, i] = np.transpose(u)
    tosave[0, i] = tk.microseconds


print(tosave)







