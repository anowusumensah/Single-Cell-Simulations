#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 17:07:21 2022

@author: cce2022
"""
import numpy as np
from scipy.integrate import solve_ivp
from AslanidiSleiman import Aslanidi
import matplotlib.pyplot as plt
import time


# This script is to run AslanidiSleiman.py script
ic=[-85.49190000, 
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

# Solving differential equations
tsart = 1000
tend = 10000
t_span = (0.0, tend)
t_eval = np.arange(0.0, tend, 0.01)
bcl = 800  #Basic cycle length
stimtimes = np.arange(tsart,tend,bcl)
mthd = 'LSODA' #LSODA is the fastest among the solvers
args = (stimtimes,)



#Solution to differential equation
start = time.time()
sol = solve_ivp(Aslanidi, t_span,ic, method = mthd, t_eval = t_eval,args = args,max_step = 0.5)
end = time.time()
print("The time elapsed", end - start)
t = sol.t
V = sol.y[0]


plt.plot(t,V)










