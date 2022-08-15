#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 11:34:10 2021
Chemostat model and Pearsons Moment Correlation
@author: cogan
"""



from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.stats import qmc # for Latin Hyptercube Sampling
from scipy.integrate import solve_ivp
import scipy as scipy
pyplot.close('all')
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible


# =============================================================================
# Define the RHS of ODEs
# =============================================================================

def f1(t,Y,params):
    N,B=Y
    N0, F, Yield, mu, K_n = params
    return N0*F-1/Yield*mu*N/(K_n+N)*B-F*N
def f2(t,Y,params):
    N,B=Y
    N0, F, Yield, mu, K_n = params
    return mu*N/(K_n+N)*B-F*B
def rhs(t,Y,params):
    N,B=Y
    N0, F, Yield, mu, K_n = params
    return [f1(t,Y,params),f2(t,Y,params)]
N_0=0
B_0=.05
# =============================================================================
# Coexistence example
# =============================================================================
params=[1, .05, .25, .5, 1]
tspan = np.linspace(0, 100, 500)
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [N_0, B_0], method='RK45',t_eval=tspan)
pyplot.figure(1)
pyplot.plot(tspan,yp.y[0])
pyplot.plot(tspan,yp.y[1])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Concentration', fontsize=16)
pyplot.legend(['Nutrient', 'Bacteria'], fontsize=16)

# =============================================================================
# Washout example
# =============================================================================
pyplot.figure(2)
params=[1, .3, .25, .5, 1]
tspan = np.linspace(0, 100, 500)
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [N_0, B_0], method='RK45',t_eval=tspan)
pyplot.figure(1)
pyplot.plot(tspan,yp.y[0])
pyplot.plot(tspan,yp.y[1])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Concentration', fontsize=16)
pyplot.legend(['Nutrient', 'Bacteria'], fontsize=16)
pyplot.savefig('Chemo_wash.png')


# =============================================================================
# Pearsons correlation coefficient
# =============================================================================
Num_samples=50
total_bacteria=np.zeros(Num_samples)

# =============================================================================
# Note that qmc is defined on the unit interval
# =============================================================================
sampler = qmc.LatinHypercube(d=5) #Define the sampling method and parameter dimension
parameters = sampler.random(n=Num_samples)      #number of samples to take
# =============================================================================
# Scale the samples into the correct parameter scale
# =============================================================================
l_bounds = np.multiply(params,.95)
u_bounds = np.multiply(params,1.05)
parameters_scaled=qmc.scale(parameters, l_bounds, u_bounds)
for i in np.arange(Num_samples):
    params_i=parameters_scaled[i,:]
    tspan = np.linspace(0, 100, 500)
    yp= solve_ivp(lambda t,Y: rhs(t,Y,params_i), [tspan[0],tspan[-1]], [N_0, B_0], method='RK45',t_eval=tspan)
# =============================================================================
# trapz is a relatively standard  implementation of the 
# trapezoidal rule for integration. QoI is total bacterial count
# =============================================================================
    total_bacteria[i]=np.trapz(yp.y[1],tspan)         
cc=np.zeros([2,5])
for j in np.arange(5):
    #Calling scipy pearsons
    cc[:,j]=scipy.stats.pearsonr(parameters_scaled[:,j], total_bacteria) 

# =============================================================================
# Increase the number of samples
# =============================================================================
Num_samples=500
total_bacteria=np.zeros(Num_samples)
sampler = qmc.LatinHypercube(d=5) 
parameters = sampler.random(n=Num_samples)      
# =============================================================================
# Scale the samples into the correct parameter scale
# =============================================================================
l_bounds = np.multiply(params,.95)
u_bounds = np.multiply(params,1.05)
parameters_scaled=qmc.scale(parameters, l_bounds, u_bounds)
for i in np.arange(Num_samples):
    params_i=parameters_scaled[i,:]
    tspan = np.linspace(0, 100, 500)
    yp_500= solve_ivp(lambda t,Y: rhs(t,Y,params_i), [tspan[0],tspan[-1]], [N_0, B_0], method='RK45',t_eval=tspan)
    total_bacteria[i]=np.trapz(yp_500.y[1],tspan)          
cc_500=np.zeros([2,5])
for j in np.arange(5):
    cc_500[:,j]=scipy.stats.pearsonr(parameters_scaled[:,j], total_bacteria) 


# =============================================================================
# Even more samples
# =============================================================================
Num_samples=1500
total_bacteria=np.zeros(Num_samples)
#Note that qmc is on the unit interval
sampler = qmc.LatinHypercube(d=5) #Define the sampling method and parameter dimension
parameters = sampler.random(n=Num_samples)      #number of samples to take
#Scale the samples into the correct parameter scale
l_bounds = np.multiply(params,.95)
u_bounds = np.multiply(params,1.05)
parameters_scaled=qmc.scale(parameters, l_bounds, u_bounds)
for i in np.arange(Num_samples):
    params_i=parameters_scaled[i,:]
    tspan = np.linspace(0, 100, 500)
    yp_1500= solve_ivp(lambda t,Y: rhs(t,Y,params_i), [tspan[0],tspan[-1]], [N_0, B_0], method='RK45',t_eval=tspan)
    total_bacteria[i]=np.trapz(yp_1500.y[1],tspan)          #trapz is a relatively standard  implementation of the trapezoidal rule for integration
cc_1500=np.zeros([2,5])
for j in np.arange(5):
    cc_1500[:,j]=scipy.stats.pearsonr(parameters_scaled[:,j], total_bacteria) 


pyplot.figure()
bar_width=.35
p_names=('$N_0$', '$F$', '$Y$', r'$\mu$', '$K_n$' )
pyplot.bar(np.arange(len(p_names))-1*bar_width,cc[0], bar_width)
pyplot.bar(np.arange(len(p_names)),cc_500[0], bar_width)
pyplot.bar(np.arange(len(p_names))+1*bar_width,cc_1500[0], bar_width)
pyplot.xticks(np.arange(len(p_names)), p_names)
pyplot.legend(['50 Samples', '500 Samples','1500 Samples'])

pyplot.savefig('cc_chemo.png')

# =============================================================================
# Don't forget to check linearity
# =============================================================================
pyplot.figure()

pyplot.scatter(parameters_scaled[:,0],total_bacteria)
pyplot.xlabel('$N_0$', fontsize=16)
pyplot.ylabel('QoI'   , fontsize=16)
pyplot.savefig('N_0_scatter.png')

pyplot.figure()

pyplot.scatter(parameters_scaled[:,1],total_bacteria)
pyplot.xlabel('$F$', fontsize=16)
pyplot.ylabel('QoI'   , fontsize=16)
pyplot.savefig('F_scatter.png')

pyplot.figure()


pyplot.scatter(parameters_scaled[:,2],total_bacteria)
pyplot.ylabel('$N_0$', fontsize=16)
pyplot.xlabel('QoI'   , fontsize=16)
pyplot.savefig('Y_scatter.png')

pyplot.figure()


pyplot.scatter(parameters_scaled[:,3],total_bacteria)
pyplot.xlabel(r'$\mu$', fontsize=16)
pyplot.ylabel('QoI'   , fontsize=16)
pyplot.savefig('muscatter.png')

pyplot.figure()


pyplot.scatter(parameters_scaled[:,4],total_bacteria)
pyplot.xlabel('$K_n$', fontsize=16)
pyplot.ylabel('QoI'   , fontsize=16)
pyplot.savefig('K_scatter.png')



