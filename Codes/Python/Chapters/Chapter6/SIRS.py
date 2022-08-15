#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 11:34:10 2021
SIRS model with cobwebbing
@author: cogan
"""



from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.integrate import solve_ivp
pyplot.close('all')
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible


#Define the RHS of ODEs
def f1(t,Y,params):
    S,I,R=Y
    k, alpha , gamma, S0, I0, R0= params
    return -k*S*I+alpha*R
def f2(t,Y,params):
    S,I,R=Y
    k, alpha , gamma, S0, I0, R0= params
    return k*S*I-gamma*I
def f3(t,Y,params):
    S,I,R=Y
    k, alpha , gamma, S0, I0, R0= params
    return gamma*I-alpha*R
def rhs(t,Y,params):
    S,I,R=Y
    k, alpha , gamma, S0, I0, R0= params
    return [f1(t,Y,params),f2(t,Y,params),f3(t,Y,params)]
p_names=( r'$k$',r'$\alpha$', r'$\gamma$', r'$S_0$',r'$I_0$',r'$R_0$' )

#With no reversion from recovered to susceptible  (alpha=0) there is only one wave
params=np.array([ 0.001, 0, .2, 500, 10, 0])
tspan = np.linspace(0, 100, 500)
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [params[-3],params[-2],params[-1]], method='RK45',t_eval=tspan)
pyplot.figure(1)

pyplot.plot(tspan,yp.y[0])
pyplot.plot(tspan,yp.y[1])
pyplot.plot(tspan,yp.y[2])
pyplot.legend(['Susceptible','Infected','Recovered'], fontsize=16)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Population', fontsize=16)
pyplot.savefig('SIRS_0.png')


#Nonzero alpha leades to oscillating solutions
params=np.array([ 0.001, 0.0125, .2, 500, 10, 0])
tspan = np.linspace(0, 200, 500)
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [params[-3],params[-2],params[-1]], method='RK45',t_eval=tspan)
pyplot.figure(2)

pyplot.plot(tspan,yp.y[0])
pyplot.plot(tspan,yp.y[1])
pyplot.plot(tspan,yp.y[2])
pyplot.legend(['Susceptible','Infected','Recovered'], fontsize=16)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Population', fontsize=16)
pyplot.savefig('SIRS_1.png')

#With a faster rate of waning antigens, the disease rapidly becomes endemic
params=np.array([ 0.001, 0.25, .2, 500, 10, 0])
tspan = np.linspace(0, 200, 500)
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [params[-3],params[-2],params[-1]], method='RK45',t_eval=tspan)
pyplot.figure(3)


pyplot.plot(tspan,yp.y[0])
pyplot.plot(tspan,yp.y[1])
pyplot.plot(tspan,yp.y[2])
pyplot.legend(['Susceptible','Infected','Recovered'], fontsize=16)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Population', fontsize=16)
#pyplot.savefig('SIRS_2.png')


# =============================================================================
# Set-up cobwebbing 
# =============================================================================

# =============================================================================
# set-up samples
# =============================================================================
params=np.array([ 0.001, 0.01, .2, 500, 10, 1])
tspan = np.linspace(0, 50, 500)
Num_samples=200
param_level=np.linspace(-1, 1, Num_samples) #this iterates the scaling of the parameter samples
parameter_list=np.zeros([Num_samples,7])
# =============================================================================
# Looop through each parameter and each value -- arrange in a matrix, use this
# to generate samples.  Note for loops are pretty fast in 
# python, but typically not in Matlab. We could make this more efficient. 
# =============================================================================
for k in np.arange(6):
    for i in np.arange(Num_samples):
        parameter_list[i,k]=(1+np.random.choice(param_level,1))*params[k]
        
for i in np.arange(Num_samples):       
    yp= solve_ivp(lambda t,Y: rhs(t,Y,parameter_list[i,0:6]), [tspan[0],tspan[-1]], [params[-3],params[-2],params[-1]], method='RK45',t_eval=tspan)
    parameter_list[i,-1]=yp.y[1][-1]
# =============================================================================
# Find all QoI (here it is the infected population at time 200) larger than a threshold
# (here 150) and sort samples into large and small infected populations
# =============================================================================
Large=np.argwhere(parameter_list[:,-1]>150) #find all QoI>
Small=np.argwhere(parameter_list[:,-1]<1) #find all QoI>
#scale parameters
parameter_list_scaled=np.zeros([Num_samples,7])
scaled=np.empty(6)

# =============================================================================
# Scale parameters between -1 and 1
# =============================================================================
for i in np.arange(Num_samples):
    np.divide(parameter_list[i,0:6], np.append(params[0:5],1), scaled)
    parameter_list_scaled[i,0:6]=scaled-1
    
# =============================================================================
#     Generate Spider plot
# =============================================================================
pyplot.figure(4)
for i in np.arange(len(Large)-1):
    pyplot.plot(parameter_list_scaled[Large[i],0:6][0],'o-',color='#FFBC79')

pyplot.plot(parameter_list_scaled[Large[-1],0:6][0],'o-',color='#FFBC79',label='High')
   

for i in np.arange(len(Small)-1):
    pyplot.plot(parameter_list_scaled[Small[i],0:6][0],'*-',color='#5F9ED1')

pyplot.plot(parameter_list_scaled[Small[-1],0:6][0],'*-',color='#5F9ED1',label='Low')
    
pyplot.xticks(ticks=np.arange(6), labels=p_names)
pyplot.legend(loc='best')
pyplot.savefig('SIRS_cobweb.png')         

