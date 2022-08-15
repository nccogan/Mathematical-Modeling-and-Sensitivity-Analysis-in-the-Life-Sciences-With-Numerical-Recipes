#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 11:34:10 2021
Noble model using Morris screening as implemented in SAlib. 
A version without packages is shown in the MATLAB implementation.
@author: cogan
"""



from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.integrate import solve_ivp
pyplot.close('all')
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible


# =============================================================================
# Define the RHS of ODEs
# =============================================================================


def f1(t,Y,params):
    V, m, h, n=Y
    Cm, alpham_0, betam_0, alphah_0, betah_0, alphan_0, betan_0,  Ek, ENa, EAn= params
    INa=(400000*m**3*h+140)*(V-ENa)
    Ik=1200*n**4*(V-Ek)+(1200*np.exp((-V-90)/50)+15*np.exp((V+90)/60))*(V-Ek)
    Ian=75*(V-ENa)
    return -(INa+Ik+Ian)/Cm
def f2(t,Y,params):
    V, m, h, n=Y
    Cm, alpham_0, betam_0, alphah_0, betah_0, alphan_0, betan_0,  Ek, ENa, EAn= params
    alpham=alpham_0*(-V-48)/(np.exp((-V-48)/15)-1)
    betam=betam_0*(V+8)/(np.exp((V+8)/5)-1)
    return alpham*(1-m)-betam*m
def f3(t,Y,params):
    V, m, h, n=Y
    Cm, alpham_0, betam_0, alphah_0, betah_0, alphan_0, betan_0,  Ek, ENa, EAn= params
    alphah=alphah_0*np.exp((-V-90)/20)
    betah=betah_0/(1+np.exp((-V-42)/10))
    return alphah*(1-h)-betah*h
def f4(t,Y,params):
    V, mh, h, n=Y
    Cm, alpham_0, betam_0, alphah_0, betah_0, alphan_0, betan_0,  Ek, ENa, EAn= params
    alphan=alphan_0*(-V-50)/(np.exp((-V-50)/10)-1)
    betan=betan_0*np.exp((-V-90)/80)
    return alphan*(1-n)-betan*n

def rhs(t,Y,params):
    V, mh, h, n=Y
    Cm, alpham_0, betam_0, alphah_0, betah_0, alphan_0, betan_0,  Ek, ENa, EAn= params
    return [f1(t,Y,params),f2(t,Y,params),f3(t,Y,params),f4(t,Y,params)]
params_list=['Cm', '$\alpha m_0$',  '$\betam_0',  '$\alphah_0',  '$\betah_0',  '$\alphan_0',  '$\betan_0',   '$\Ek',  '$\ENa',  '$\EAn']

params=[6.000, 100, 120,170,1000,0.1,2, -100, 40,-60]
V0=-85
m0=.5
h0=.75
n0=.65
tspan = np.linspace(0, 1, 500)
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [V0,m0,h0,n0], method='LSODA',t_eval=tspan)
pyplot.plot(tspan,yp.y[0])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Voltage', fontsize=16)
pyplot.savefig('noble_trace.png')

# =============================================================================
# plot histogram example
# =============================================================================
bin_number,frequency,_=pyplot.hist(np.abs((np.fft.fft(yp.y[0]))),500) 
pyplot.axis([0, 1000, 0, 100])
pyplot.xlabel('Frequency', fontsize=16)
pyplot.ylabel('Occurance Frequency', fontsize=16)
pyplot.savefig('noble_FFT.png')

# =============================================================================
# Measure the frequency using fft
# =============================================================================
bin_number,frequency = np.histogram((np.fft.fft(yp.y[0])),500)
max_bin=np.argmax(bin_number) #find bin with max frequency
QoI=np.abs(frequency[max_bin])
# =============================================================================
# Morris method as implemented in SAlib
# =============================================================================

from SALib.analyze import morris
from SALib.sample.morris import sample
from SALib.plotting.morris import horizontal_bar_plot, covariance_plot, \
    sample_histograms

params=[6.000, 100, 120,170,1000,0.1,2, -100, 40,-60]
V0=-85
m0=.5
h0=.75
n0=.65
b=np.empty([10,2])

# =============================================================================
# Problem defines the sampling space. It is a python library that defines 
# variables in ''. 
# =============================================================================
for i in np.arange(10):
    b[i,:]=[.95*params[i],1.05*params[i]]
problem = {
  'num_vars': 10,
  'names': [r'Cm', r'${\alpha_m}_0$',  r'${\beta_m}_0$',  r'${\alpha_h}_0$',  r'${\beta_h}_0$',  r'${\alpha_n}_0$',  r'${\beta_n}_0$',   r'$E_k$',  r'$E_{Na}$',  r'$E_{An}$'],
  'groups': None,
  'bounds': [[ 5.700e+00,  6.300e+00],
         [ 9.500e+01,  1.050e+02],
         [ 1.140e+02,  1.260e+02],
         [ 1.615e+02,  1.785e+02],
         [ 9.500e+02,  1.050e+03],
         [ 9.500e-02,  1.050e-01],
         [ 1.900e+00,  2.100e+00],
         [-1.050e+02, 9.500e+01],
         [ 3.800e+01,  4.200e+01],
         [ -6.300e+01,-5.700e+01]]
}

# =============================================================================
#  Generate samples using the sample command
# =============================================================================
param_values = sample(problem, N=500, num_levels=4,
                      optimal_trajectories=None)

# =============================================================================
#  Run the "model" -- this will happen offline for external models
# =============================================================================
tspan = np.linspace(0, 1, 10)
Y = np.empty([len(param_values)])
for i in np.arange(len(Y)):
    yp= solve_ivp(lambda t,Y: rhs(t,Y,param_values[i]), [tspan[0],tspan[-1]], [V0,m0,h0,n0], method='LSODA',t_eval=tspan)
    bin_number,frequency = np.histogram(np.abs(np.fft.fft(yp.y[0])),500)
    max_bin=np.argmax(bin_number) #find bin with max frequency
    Y[i]=np.abs(frequency[max_bin])
 
# =============================================================================
#  Perform the sensitivity analysis using the model output
#  Specify which column of the output file to analyze (zero-indexed)
# =============================================================================
Si = morris.analyze(problem, param_values, Y, conf_level=0.95,
                    print_to_console=True,
                    num_levels=4, num_resamples=100)
# =============================================================================
#  Returns a dictionary with keys 'mu', 'mu_star', 'sigma', and 'mu_star_conf'
#  e.g. Si['mu_star'] contains the mu* value for each parameter, in the
#  same order as the parameter file
# =============================================================================
fig, (ax1) = pyplot.subplots(1, 1)
horizontal_bar_plot(ax1,Si, {}, sortby='mu_star', unit=r"QoI")
#covariance_plot(ax2, Si, {}, unit=r"tCO$_2$/year")

