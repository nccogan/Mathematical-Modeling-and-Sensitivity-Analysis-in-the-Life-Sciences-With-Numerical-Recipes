#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 07:30:47 2022
Solves the SI model with vital dynamics. Uses Spider plot
as the SI indice.
@author: cogan
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
from matplotlib import pyplot    # import plotting library
from scipy.integrate import solve_ivp
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible

#Define the RHS of ODEs
def f1(t,Y,params):
    S,I=Y
    r, kappa, k, delta, IC1, IC2 = params
    return -k*S*I+r*S*(kappa-S)
def f2(t,Y,params):
    S,I=Y
    r, kappa, k, delta , IC1, IC2= params
    return k*S*I-delta*I
def rhs(t,Y,params):
    S,I=Y
    r, kappa, k, delta , IC1, IC2= params
    return [f1(t,Y,params),f2(t,Y,params)]
params_names=( ['$r$', '$\kappa$', '$k$', '$\delta$', '$S_0$', '$I_0$'])

params=np.array([.05, 100, .075, .3,99,1])
tspan = np.linspace(0, 25, 500)
y_solution= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [99,1], method='RK45',t_eval=tspan)
pyplot.plot(tspan, y_solution.y[0])
pyplot.plot(tspan, y_solution.y[1])
pyplot.legend(['Susceptible', 'Infected'], fontsize=16)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Population', fontsize=16)
#pyplot.savefig('SI_dyn.png')

#####################################
tspan = np.linspace(0, 500, 1000)

params=[.05, 100, .075, .3,99,1] #augment parameters for initial conditions
Num_samples=5
QoI=np.zeros([2*Num_samples,5])
for k in np.arange(0,5): 
    for i in np.arange(0,2*Num_samples):
        params[k]=.5*params[k]+i*.1*params[k]
        y0=[params[-2],params[-1]]
        y_solution = solve_ivp(lambda t,Y: rhs(t,Y,params),  [tspan[0],tspan[-1]],y0,method='RK45',t_eval=tspan)
        QoI[i,k]=y_solution.y[0][-1]/(y_solution.y[0][-1]+y_solution.y[1][-1])
        params=[.05, 100, .075, .3,99,1]
pyplot.figure()
pyplot.plot(QoI,'-o')
pyplot.legend(params_names, fontsize=16)
pyplot.xticks(np.arange(0,2*Num_samples+1), ['-50%', '-40%', '-30%', '-20%', '-10%', '0%', '10%', '20%', '30%', '40%', '50%'])
#pyplot.savefig('SI_spider.png')
