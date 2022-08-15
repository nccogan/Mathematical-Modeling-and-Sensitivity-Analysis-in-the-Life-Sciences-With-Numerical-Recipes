#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 11:34:10 2021
FHN model 
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
    V, w=Y
    epsilon, A, alpha, w0, gamma= params
    return 1/epsilon*(A*V*(V-alpha)*(1-V)-w-w0) 
def f2(t,Y,params):
    V, w=Y
    epsilon, A, alpha, w0, gamma= params
    return (V-gamma*w)

def rhs(t,Y,params):
    V, w=Y
    epsilon, A, alpha, w0, gamma= params
    return [f1(t,Y,params),f2(t,Y,params)]
params_list=[r'$\epsilon$', '$A$', r'$\alpha$', '$w_0$', r'$\gamma$']
#Change IC to show excitable
params=[.01, 1, .1, 0, .5]
V0=.1
w0=0
tspan = np.linspace(0, 2, 500)
yp1= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [.1,w0], method='LSODA',t_eval=tspan)
yp2= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [.2,w0], method='LSODA',t_eval=tspan)
pyplot.plot(tspan,yp1.y[0])
pyplot.plot(tspan,yp2.y[0])
pyplot.legend(['$V_0=.1$', '$V_0=.2$'],fontsize=16)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Voltage', fontsize=16)
pyplot.savefig('FHN_Excitable.png')


params=[.01, 1, .1, -1, .5]
V0=.1
w0=0
tspan = np.linspace(0, 10, 500)
yp1= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [.1,w0], method='LSODA',t_eval=tspan)
pyplot.plot(tspan,yp1.y[0])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Voltage', fontsize=16)
pyplot.savefig('FHN_tonic.png')

# =============================================================================
# sensitivity section
# =============================================================================
# =============================================================================
# First set-up simulate data
# =============================================================================
params=[.01, 1, .1, -1, .5]
V0=.1
w0=0
tspan = np.linspace(0, 10, 50)
y_soln= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [V0,w0], method='LSODA',t_eval=tspan)
y_data=y_soln.y[0]+.1*np.random.uniform(-1.0, 1.0,len(y_soln.y[0]))
pyplot.plot(y_data,'o')
pyplot.plot(y_soln.y[0],'*')
pyplot.legend(['True', 'With Noise'])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Voltage', fontsize=16)


# =============================================================================
# Begin Sampling and imports
# =============================================================================
from SALib.analyze import sobol
from SALib.analyze import delta

from SALib.sample.morris import sample
#from SALib.plotting.morris import horizontal_bar_plot, covariance_plot, \
#    sample_histograms

# or define parameters for sampling
b=np.empty([len(params),2])

for i in np.arange(len(params)):
    b[i,:]=[.95*params[i],1.05*params[i]]
b.sort(axis=1) # sort the rows that are negative    

problem = {
  'num_vars': 5,
  'names': [r'$\epsilon$', '$A$', r'$\alpha$', '$w_0$', r'$\gamma$'],
  'groups': None,
  'bounds': b
}
# Files with a 4th column for "group name" will be detected automatically, e.g.
# param_file = '../../src/SALib/test_functions/params/Ishigami_groups.txt'

# Generate samples
param_values = sample(problem, N=1000, num_levels=4,
                      optimal_trajectories=None)

# To use optimized trajectories (brute force method),
# give an integer value for optimal_trajectories

# Run the "model" -- this will happen offline for external models
Y = np.empty([len(param_values)])
for i in np.arange(len(Y)):
    yp= y_soln= solve_ivp(lambda t,Y: rhs(t,Y,param_values[i]), [tspan[0],tspan[-1]], [V0,w0], method='LSODA',t_eval=tspan)
    Y[i]=np.linalg.norm(yp.y[0]-y_data)
# Perform the sensitivity analysis using the model output
# Specify which column of the output file to analyze (zero-indexed)
Si = sobol.analyze(problem, Y, calc_second_order=True, conf_level=0.95, print_to_console=True)
Si1 = delta.analyze(problem, param_values, Y, print_to_console=True)
bar_width=.35
pyplot.bar(np.arange(len(params_list))-.5*bar_width,Si['S1'], bar_width)
pyplot.bar(np.arange(len(params_list))+.5*bar_width,Si1['S1'], bar_width)
pyplot.legend(['Sobol', 'Moment Independent'],fontsize=16)
pyplot.xticks(np.arange(len(params_list)), params_list)
pyplot.ylabel('SA',fontsize=16)

pyplot.savefig('FHN_delta.png')
# =============================================================================
# compare with varying parameters to compare variations in output. 
# One could be more precise by sampling over reduced parameter space.
# =============================================================================
params=[.01, 1, .1, -1, .5]
V0=.1
w0=0
tspan = np.linspace(0, 10, 500)
y_true= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [V0,w0], method='LSODA',t_eval=tspan)
params=[.01, 1, 1, -1, .5]
V0=.1
w0=0
tspan = np.linspace(0, 10, 500)
y_alpha= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [V0,w0], method='LSODA',t_eval=tspan)
params=[.1, 1, .1, -1, .5]
V0=.1
w0=0
tspan = np.linspace(0, 10, 500)
y_eps= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [V0,w0], method='LSODA',t_eval=tspan)
params=[.01, 1, .1, -1.5, .5]
V0=.1
w0=0
tspan = np.linspace(0, 10, 500)
y_w= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [V0,w0], method='LSODA',t_eval=tspan)
pyplot.plot(tspan, y_true.y[0])
pyplot.plot(tspan, y_alpha.y[0])
pyplot.legend(['True', r'$\alpha=.5$'],fontsize=16)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Voltage', fontsize=16)
pyplot.savefig('FHN_compare_alpha.png')


pyplot.plot(tspan, y_true.y[0])
pyplot.plot(tspan, y_eps.y[0])
pyplot.legend(['True', r'$\epsilon=.1$'],fontsize=16)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Voltage', fontsize=16)
pyplot.savefig('FHN_compare_eps.png')


pyplot.plot(tspan, y_true.y[0])
pyplot.plot(tspan, y_w.y[0])
pyplot.legend(['True', r'$w_0=-1.5$'],fontsize=16)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Voltage', fontsize=16)
pyplot.savefig('FHN_compare_w.png')
