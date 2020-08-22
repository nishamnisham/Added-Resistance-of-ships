# -*- coding: utf-8 -*-
"""
@author: nisham
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import math

#ship specifications
lbp=230
b=32
t=10.8
cb=0.61
le=64.778
ro=1025
zeta=1
kyy = 0.25*lbp
amp=1
dtr=0

lamdaratio = np.arange(0.2,3.1,0.05)
lmda = np.array([num*lbp for num in lamdaratio])


def rawref(v):
    '''Added wave resistance due to wave refraction is calculated
    throughout this function. The relevant variables are found within'''
    fn= (v*.51444)/(math.sqrt(9.81*lbp))
    
    E = math.atan(b/(2*le))
    
    alphaT=[]
    for i in range(len(lamdaratio)):
        if (lamdaratio[i] <= 2.5).any():
            alphaT.append(1-(math.exp(((-4*math.pi)*\
                                       ((t/lmda[i])-(t/(2.5*lbp)))))))
        elif (lamdaratio[i] > 2.5).any():
            alphaT.append(0)
    alphaT=np.array(alphaT)
            
    rawr =[]
    for i in range(len(lamdaratio)):
        rawr.append(((2.25/2)*ro*9.81*b*(amp**2)*(alphaT[i]))*\
                    (math.sin(E)**2)*(1+(5*math.sqrt(lbp/lmda[i])*fn))*\
                    ((0.87/cb)**(1+(4.0*(math.sqrt(fn))))))
    rawr = np.array(rawr)

    return rawr


def rawmotion(v):
    '''Added resistance due to motion induced by regular waves
    are calculated within this function. Note some global parameters 
    are used here. This may not work without providing them first.'''
    fn= (v*.51444)/(math.sqrt(9.81*lbp))
    
    #a1
    a1 = (60.3*(cb**1.34))*((4*kyy/lbp)**2)*\
    ((0.87/cb)**(1+fn))*((math.log(b/t))**(-1))

    #alpha2
    a2 = (fn**1.5)*(math.exp(-3.5*fn))
    #alpha3
    a3 = 1.0+(0.25*(math.atan(dtr/lbp)))
    #Omega Calculations
    omega =[]
    for i in range(len(lamdaratio)):
        omega.append((2.142*((kyy/lbp)**0.333))*(math.sqrt(lbp/lmda[i]))*\
                     (1-(0.13*(0.85/cb)*\
                         (math.log(b/t)-math.log(2.75))*(fn**0.143))))
    omega=np.array(omega)
    b1 =[]
    d1=[]
    if cb <= 0.75:
        for i in range(len(omega)):
            if (omega[i]<1).any():
                b1.append(11.0)
                d1.append(14)
            else:
                b1.append(-8.5)
                d1.append((-566*((lbp/b)**(-2.66)))*6)
    else:
        for i in range(len(omega)):
            if omega[i]<1:
                b1.append(11.0)
                d1.append(-566*(lbp/b)**(-2.66))
            else:
                b1.append(-8.5)
                d1.append(-566*((lbp/b)**(-2.66))*6)
    b1 = np.array(b1)
    d1 = np.array(d1)     
    rawm =[]
    for i in range(len(lamdaratio)):
        rawm.append((4*ro*9.81*(amp**2)*(b**2)*(omega[i]**(b1[i])))*\
                    math.exp(((b1[i]/d1[i])*(1-omega[i]**(d1[i]))))*\
                    (a1*a2*a3/lbp))
    rawm =np.array(rawm)
    return rawm

#desired ship speeds
vs = np.array([18,24])

rawrs = np.array([rawref(v) for v in vs])
rawms = np.array([rawmotion(v) for v in vs])
raw = np.add(rawrs,rawms)

def Rdg():
    Rdgv1 = np.array([num/(ro*9.81*(amp**2)*(b**2)/lbp) for num in raw[0,0:]])
    Rdgv2 = np.array([num/(ro*9.81*(amp**2)*(b**2)/lbp) for num in raw[1,0:]])
    return Rdgv1,Rdgv2

x = lamdaratio
y1= Rdg()[0]
y2= Rdg()[1]

line1 = plt.plot(x,y1, color='blue')
line1 = plt.plot(x,y2, color='red')
plt.xlabel('lamda/L')
plt.ylabel('Raw/ro.g.a2.b2/lbp')
plt.show()