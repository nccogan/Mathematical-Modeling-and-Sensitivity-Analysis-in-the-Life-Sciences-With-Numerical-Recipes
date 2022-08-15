#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 11:34:10 2021
Chemostat model : Steady-state
@author: cogan
"""


#### 
from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
pyplot.close('all')
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible


# =============================================================================
# Define the RHS of ODEs for nutrient (N) and bacteria (B)
# =============================================================================
def f1(Y,params):
    N,B=Y
    N0, F, Yield, mu, K_n = params
    return N0*F-1/Yield*mu*N/(K_n+N)*B-F*N
def f2(Y,params):
    N,B=Y
    N0, F, Yield, mu, K_n = params
    return mu*N/(K_n+N)*B-F*B
def rhs(Y,params):
    N,B=Y
    N0, F, Yield, mu, K_n = params
    return [f1(Y,params),f2(Y,params)]
N_0=1
B_0=.05
params=[1, .05, .25, .5, 1]
fp = []
solution1=fsolve(lambda Y: rhs(Y,params), [.1, .1])
print('Fixed points = ',solution1)
# =============================================================================
# A discrete estimate of the Jacobian using centered difference. The user should
#consider this with care.
# =============================================================================
def Myjac(Y,params,epsilon):
    N,B=Y
    f_N = (f1([N+epsilon,B],params)-f1([N-epsilon,B],params))/epsilon/2
    f_B = (f1([N,B+epsilon],params)-f1([N,B-epsilon],params))/epsilon/2
    g_N = (f2([N+epsilon,B],params)-f2([N-epsilon,B],params))/epsilon/2
    g_B = (f2([N,B+epsilon],params)-f2([N,B-epsilon],params))/epsilon/2
    J = np.array([[f_N,f_B],[g_N,g_B]])
    return J

J=Myjac(np.array([solution1[0],solution1[1]]),params,.001)
print('Jacobian is', J)

# =============================================================================
# Determine the eigenvalues
# =============================================================================
w,v= np.linalg.eig(J)
print('Eigenvalues are', w)
tspan = np.linspace(0, 2000, 4000)

y_soln= solve_ivp(lambda t,Y: rhs(Y,params), [tspan[0],tspan[-1]], [1.10,.1], method='LSODA',t_eval=tspan)
pyplot.plot(tspan,y_soln.y[0])
pyplot.plot(tspan,y_soln.y[1])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Concentration', fontsize=16)
pyplot.legend(['Nutrient', 'Bacteria'], fontsize=16)
