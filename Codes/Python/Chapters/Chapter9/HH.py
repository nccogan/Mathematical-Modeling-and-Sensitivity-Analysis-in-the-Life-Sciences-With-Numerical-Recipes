#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 11:34:10 2021
Hodgkin-Huxley model: Sobol Sensitivity as implemented in SAlib
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
    V, m, h, n=Y
    Cm, alpham_0, betam_0, alphah_0, betah_0, alphan_0, betan_0,gK,gNa, gl, Vk, VNa, Vl, Iapp= params
    INa=gNa*m**3*h*(V-VNa)
    Ik=gK*n**4*(V-Vk)
    Il=gl*(V-Vl)
    return -(INa+Ik+Il)/Cm+Iapp
def f2(t,Y,params):
    V, m, h, n=Y
    Cm, alpham_0, betam_0, alphah_0, betah_0, alphan_0, betan_0,gK,GNa, Gl, Vk, VNa, Vl, Iapp= params
    alpham=alpham_0*(25-V)/(np.exp((25-V)/10)-1)
    betam=betam_0*np.exp(-V/18)
    return alpham*(1-m)-betam*m
def f3(t,Y,params):
    V, m, h, n=Y
    Cm, alpham_0, betam_0, alphah_0, betah_0, alphan_0, betan_0,gK,GNa, Gl, Vk, VNa, Vl, Iapp= params
    alphah=alphah_0*np.exp((-V)/20)
    betah=betah_0/(1+np.exp((-V+30)/10))
    return alphah*(1-h)-betah*h
def f4(t,Y,params):
    V, mh, h, n=Y
    Cm, alpham_0, betam_0, alphah_0, betah_0, alphan_0, betan_0,gK,GNa, Gl, Vk, VNa, Vl, Iapp= params
    alphan=alphan_0*(-V+10)/(np.exp((-V+10)/10)-1)
    betan=betan_0*np.exp((-V)/80)
    return alphan*(1-n)-betan*n

def rhs(t,Y,params):
    V, mh, h, n=Y
    Cm, alpham_0, betam_0, alphah_0, betah_0, alphan_0, betan_0,gK,GNa, Gl, Vk, VNa, Vl, Iapp= params
    return [f1(t,Y,params),f2(t,Y,params),f3(t,Y,params),f4(t,Y,params)]
params_list=[ r'$C_m$', r'${\alpha_m}_0$', r'${\beta_m}_0$', r'${\alpha_h}_0$', r'${\beta_h}_0$', r'${\alpha_n}_0$', r'${\beta_n}_0$',r'$g_K$',r'$g_{Na}$', r'$g_l$', r'$V_k$', r'$V_{Na}$', r'$V_l$', r'$I_{app}$']


params=[2.000, .1, 4,.07,1,0.01,.125,36, 120, .3,-12, 115,10.6,0*35]
V0=1
m0=.053
h0=.595
n0=.317
tspan = np.linspace(0, 20, 500)
yp1= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [V0,m0,h0,n0], method='LSODA',t_eval=tspan)
params=[2.000, .1, 4,.07,1,0.01,.125,36, 120, .3,-12, 115,10.6,0*35]
V0=6
m0=.053
h0=.595
n0=.317
tspan = np.linspace(0, 20, 500)
yp2= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [V0,m0,h0,n0], method='LSODA',t_eval=tspan)
params=[2.000, .1, 4,.07,1,0.01,.125,36, 120, .3,-12, 115,10.6,0*35]
V0=20
m0=.053
h0=.595
n0=.317
tspan = np.linspace(0, 20, 500)
yp3= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [V0,m0,h0,n0], method='LSODA',t_eval=tspan)

pyplot.plot(tspan,yp1.y[0])
pyplot.plot(tspan,yp2.y[0])

pyplot.plot(tspan,yp3.y[0])
pyplot.legend(['$V_0=1$','$V_0=6$','$V_0=20$'])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Voltage', fontsize=16)
pyplot.savefig('HH_compare.png')
# =============================================================================
# Sobol' sensitivity. The setup is essentially the same as Morris screening
# since they are both calculated using the the SAlib package
# =============================================================================

from SALib.analyze import sobol
from SALib.sample.morris import sample

params=[2.000, .1, 4,.07,1,0.01,.125,36, 120, .3,-12, 115,10.6,10]
V0=1
m0=.053
h0=.595
n0=.317
# or define parameters
b=np.zeros([14,2])
# =============================================================================
# Set up the samples
# =============================================================================

for i in np.arange(len(params)):
    b[i,:]=[.95*params[i],1.05*params[i]]
b.sort(axis=1) # sort the rows that are negative    
problem = {
  'num_vars': 13,
  'names': [ r'$C_m$', r'${\alpha_m}_0$', r'${\beta_m}_0$', r'${\alpha_h}_0$', r'${\beta_h}_0$', r'${\alpha_n}_0$', r'${\beta_n}_0$',r'$g_K$',r'$g_{Na}$', r'$g_l$', r'$V_k$', r'$V_{Na}$', r'$V_l$', r'$I_{app}$'],
  'groups': None,
  'bounds': b
}
# Files with a 4th column for "group name" will be detected automatically, e.g.
# param_file = '../../src/SALib/test_functions/params/Ishigami_groups.txt'

# Generate samples
param_values = sample(problem, N=900, num_levels=4,
                      optimal_trajectories=None)

# =============================================================================
# Run the model
# =============================================================================

tspan = np.linspace(0, 20, 100)
Y = np.empty([len(param_values)])
for i in np.arange(len(Y)):
    yp= solve_ivp(lambda t,Y: rhs(t,Y,param_values[i]), [tspan[0],tspan[-1]], [V0,m0,h0,n0], method='LSODA',t_eval=tspan)
    Y[i]=np.max(yp.y[0])
# =============================================================================
# Perform the sensitivity analysis using the model output
# =============================================================================
Si = sobol.analyze(problem, Y, calc_second_order=True, conf_level=0.95, print_to_console=True)
# Returns a dictionary with keys 'S1', 'S1_conf', 'ST', and 'ST_conf'
# e.g. Si['S1'] contains the first-order index for each parameter,
# in the same order as the parameter file
bar_width=.35
pyplot.bar(np.arange(len(params_list))-.5*bar_width,Si['S1'], bar_width)
pyplot.xticks(np.arange(len(params_list)), params_list)
pyplot.ylabel('Sobol',fontsize=16)
