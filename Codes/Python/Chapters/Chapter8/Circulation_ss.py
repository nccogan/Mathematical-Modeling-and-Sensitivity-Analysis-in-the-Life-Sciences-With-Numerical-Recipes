#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 11:34:10 2021
Circulation model: Algebraic 
@author: cogan
"""



from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.stats import qmc # for LHS
from scipy.optimize import fsolve
pyplot.close('all')
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible


#Define the RHS of the equations
def f1(Y,params):
    Q,Pa, Pv=Y
    Qvessel, R, gamma, Vvessel,V0, F, Vmax, Vmin, Cd, Cs= params
    return (Pa-Pv)/R*(1+gamma*(Pa+Pv)+gamma**2/3*(Pa**2+Pa*Pv+Pv**2))-Qvessel
def f2(Y,params):
    Q,Pa, Pv=Y
    Qvessel, R, gamma, Vvessel,V0, F, Vmax, Vmin, Cd, Cs= params
    return V0*(1+gamma/2*(Pa+Pv)+gamma**2/6*(Pa-Pv)**2)-Vvessel
def f3(Y,params):
    Q,Pa, Pv=Y
    Qvessel, R, gamma, Vvessel,V0, F, Vmax, Vmin, Cd, Cs= params
    return F*(Vmax-Vmin+Cd*Pv-Cs*Pa)-Q
def rhs(Y,params):
    Q,Pa, Pv=Y
    Qvessel, R, gamma, Vvessel,V0, F, Vmax, Vmin, Cd, Cs= params
    return [f1(Y,params),f2(Y,params),f3(Y,params)]

params=[5, 5.5, 1.50, 11,2,80, 5,4, .5,.12]
solution1 = fsolve(rhs, (10,.11,.1),args=(params,) )
print(solution1)

# =============================================================================
# Normal distribution
# =============================================================================
n_power=10
Num_samples=2**n_power
total_bacteria=np.zeros(Num_samples)

import pyDOE as doe
from scipy.stats.distributions import norm
parameters_normal=doe.lhs(10,samples=Num_samples)

for i in np.arange(10):
     parameters_normal[:,i] = norm(loc=params[i], scale=.04*params[i]).ppf(parameters_normal[:, i]) #scaled parameters
#Explore how the variance changes things

# =============================================================================
# Uniform distribution
# =============================================================================
sampler = qmc.LatinHypercube(d=10) #Define the sampling method and parameter dimension
parameters_uniform = sampler.random(n=Num_samples)      #number of samples to take
 #Scale the samples into the correct parameter scale
l_bounds = np.multiply(params,.95)
u_bounds = np.multiply(params,1.05)
parameters_uniform=qmc.scale(parameters_uniform, l_bounds, u_bounds)
# 
# =============================================================================
# Sobol distribution: Scrambled Sobol Sequences.
# =============================================================================
sampler = qmc.Sobol(d=10,scramble=True,) #Define the sampling method and parameter dimension
parameters_sobol = sampler.random_base2(m=n_power)   #number of samples to take
 #Scale the samples into the correct parameter scale
l_bounds = np.multiply(params,.95)
u_bounds = np.multiply(params,1.05)
parameters_sobol=qmc.scale(parameters_sobol, l_bounds, u_bounds)

FP_uniform=np.zeros([3,Num_samples])
FP_normal=np.zeros([3,Num_samples])
FP_sobol=np.zeros([3,Num_samples])


for i in np.arange(Num_samples):
    FP_uniform[:,i]=fsolve(rhs, (10,.11,.1),args=(parameters_uniform[i,:],) )
    FP_normal[:,i]=fsolve(rhs, (10,.11,.1),args=(parameters_normal[i,:],) )
    FP_sobol[:,i]=fsolve(rhs, (10,.11,.1),args=(parameters_sobol[i,:],) )


pyplot.hist(FP_uniform[0,:],20, alpha=.5)
pyplot.hist(FP_normal[0,:],20, alpha=.5)
pyplot.hist(FP_sobol[0,:],20, alpha=.5)
pyplot.axis([0, 175,0, 300])
pyplot.xlabel('$QoI=Q$', fontsize=16)
pyplot.ylabel('Frequency', fontsize=16)
pyplot.legend(['Uniform', 'Normal', 'Sobol'])
pyplot.savefig('Q.png')

pyplot.figure()
pyplot.hist(FP_uniform[1,:],20, alpha=.5)
pyplot.hist(FP_normal[1,:],20, alpha=.5)
pyplot.hist(FP_sobol[1,:],20, alpha=.5)
pyplot.axis('tight')
pyplot.xlabel('$QoI=P_a$', fontsize=16)
pyplot.ylabel('Frequency', fontsize=16)
pyplot.legend(['Uniform', 'Normal', 'Sobol'])
pyplot.savefig('P_a.png')

pyplot.figure()

pyplot.hist(FP_uniform[2,:],20, alpha=.5)
pyplot.hist(FP_normal[2,:],20, alpha=.5)
pyplot.hist(FP_sobol[2,:],20, alpha=.5)
pyplot.axis('tight')
pyplot.xlabel('$QoI=P_v$', fontsize=16)
pyplot.ylabel('Frequency', fontsize=16)
pyplot.legend(['Uniform', 'Normal', 'Sobol'])
pyplot.savefig('P_v.png')


   