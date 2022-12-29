#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 13:33:29 2022

@author: cce2022
"""
import numpy as np
import matplotlib.pyplot as plt
from VoltageClamp import VoltageProtocol
from VoltageClamp import AslanidiVoltageClamp
import time
from scipy.integrate import solve_ivp



ic=[ 
    0.00142380, 
    0.98646900,
    0.99907100,
    0.00142380,
    0.93722800,
    0.00195127,
    0.00036627,
    0.99816200,
    0.99816200,
    0.03065770,
    1.71221e-06,
    0.99764600,
    8.06835e-05,
    0.97941500,
    9.93905e-02,
    1.28293,
    0.000131111,
    0.0121366,
    0.0116031,
    0.132481,
    0.0981466,
    0.633156,
    0.00026248,
    0.04,
    2.67486e-07,
    1.94911e-06,
    0.879232,
    0.000336009,
    0.00244706,
    0.00243042,
    0.00996049,
    0.12297,
    0.136961,
    0.0079682,
    0.740619326616887
    
    ]; 

ms2s =1000
dt = 0.01
tend = 5000
t_span = (0.0, tend)
t_eval = np.arange(0.0, tend, dt)
splice = np.where(t_eval >= 2*ms2s)[0][0] #Index where depolarization voltage starts
mthd = 'BDF' #LSODA is the fastest among the solvers

#Parameters for IKr
T = 310.0; 
F     = 96486.7;        #units(J / kmol / K);
R     = 8314.3;         #units(C / mol); 
RTonF = R * T / F;          #units(mV);
GKr   = 0.015583333333333;       #units(mS / cm^2);
Ki      = 135.0;                 #units(mM);
Ke      = 5.4;                   #units(mM);
E_K  = RTonF * np.log(Ke  / Ki )  

Vcm = np.arange(40,-151,-10) #Command voltages
lenVm = Vcm.shape[0] #length of Vcm
lenTeval = t_eval.shape[0] #length of t_eval
storeIKrVcm = np.zeros([lenVm,lenTeval]) #allocate memory for IKr
storeXrVcm = np.zeros([lenVm,lenTeval]) #allocate memory for the activation parameter
for iterVc in np.arange(lenVm):
    start = time.time()
    sol = solve_ivp(AslanidiVoltageClamp, t_span,ic, method = mthd, t_eval = t_eval,args = (Vcm[iterVc],),max_step = 0.5)
    end = time.time()
    print("The time elapsed is {} for {}mV".format(end - start, Vcm[iterVc]) )
    #solve the current
    t = sol.t
    xr = sol.y[5]
    V = VoltageProtocol(t, Vcm[iterVc])
    #x_inf = 1.1/(1.0+np.exp(-1.0*(V-1.5)/20));
    r_infty = 1.0/(1.0 + np.exp(V / 50));
    IKr = GKr * xr * r_infty * (V - E_K);
    #Keep track of Xr and IKr current for the for the voltage steps
    storeXrVcm[iterVc,:] = xr 
    storeIKrVcm[iterVc,:] = IKr 
    

peakXr = np.amax(storeXrVcm[:,splice:], axis=1)
peakIKr = np.amax(storeIKrVcm[:,splice:], axis=1)    

#open probability plot
fig, ax = plt.subplots()
#plt.subplot(211)
plt.title('Open Probabilty of IKr current')
plt.plot(Vcm, peakXr)
plt.xlabel('Voltage (mV)')
plt.ylabel('Open Probability')
ax.set_xlim(Vcm[-1], Vcm[0])
plt.show()

#I-V plot
#fig, ax = plt.subplot(212)
plt.title('I-V plot')
plt.plot(Vcm, peakIKr)
plt.xlabel('Voltage (mV)')
plt.ylabel('Current (pA/pF)')
ax.set_xlim(Vcm[-1], Vcm[0])
plt.show()




for i in np.arange(lenVm):
    plt.plot(sol.t, storeIKrVcm[i,:])










