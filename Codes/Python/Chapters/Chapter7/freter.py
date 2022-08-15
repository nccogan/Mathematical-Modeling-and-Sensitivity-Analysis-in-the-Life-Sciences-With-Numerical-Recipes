
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 11:34:10 2021
Freter model 
We are using pyDOE to do some of the sampling. 
This may require the command
'pip install --upgrade pyDOE'
@author: cogan
"""



from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.integrate import solve_ivp
import scipy as scipy
pyplot.close('all')
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible


# =============================================================================
# Define the RHS of ODEs
# =============================================================================
def f1(t,Y,params):
    N, Bu, Bb=Y
    N0, Yield, mu, K_n, F, B0, alpha, K_alpha, V, A, beta, K_b = params
    return N0*F-1/Yield*mu*N/(K_n+N)*(Bu+Bb)-F*N
def f2(t,Y,params):
    N, Bu, Bb=Y
    N0, Yield, mu, K_n, F, B0, alpha, K_alpha, V, A, beta, K_b = params
    return B0*F-alpha/(K_alpha+Bb)*Bu+V/A*beta*Bb \
            -F*Bb +(1-Bb/(K_b+Bb))*1/Yield*mu*N/(K_n+N)*Bb-F*Bu
def f3(t,Y,params):
    N, Bu, Bb=Y
    N0, Yield, mu, K_n, F, B0, alpha, K_alpha, V, A, beta, K_b = params
    return A/V*alpha/(K_alpha+Bb)*Bu-beta*Bb \
        +A/V*Bb/(K_b+Bb)*1/Yield*mu*N/(K_n+N)*Bb
def rhs(t,Y,params):
    N, Bu, Bb=Y
    N0, Yield, mu, K_n, F, B0, alpha, K_alpha, V, A, beta, K_b = params
    return [f1(t,Y,params),f2(t,Y,params),f3(t,Y,params)]
Params_names=['$N_0$', '$Y$', '$\mu$', '$K_n$', '$F$', '$B_0$', r'$\alpha$', r'$k_\alpha$', '$V$', '$A$', r'$\beta$', '$K_b$'] 
# =============================================================================
# Coexistence
# =============================================================================
params=[1, .1, 1, .5, 1, 0*.1, .1, .1, 10, 1, .2, .5]
tspan = np.linspace(0, 40, 500)
N_init=0
Bu_init=.1
Bb_init=0
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [N_init, Bu_init, Bb_init], method='RK45',t_eval=tspan)
pyplot.figure(1)
pyplot.plot(tspan,yp.y[0])
pyplot.plot(tspan,yp.y[1])
pyplot.plot(tspan,yp.y[2])

pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Concentration', fontsize=16)
pyplot.legend(['$N$', '$B_u$', '$B_b$'], fontsize=16)
pyplot.savefig('Freter_coexist.png')

pyplot.figure(2)

# =============================================================================
# Washout
# =============================================================================
params=[1, .1, .1, .5, 2, 0*.1, .1, .1, 10, 1, .2, .5]
tspan = np.linspace(0, 4, 500)
N_init=0
Bu_init=.1
Bb_init=0
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [N_init, Bu_init, Bb_init], method='RK45',t_eval=tspan)
pyplot.figure(1)
pyplot.plot(tspan,yp.y[0])
pyplot.plot(tspan,yp.y[1])
pyplot.plot(tspan,yp.y[2])

pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Concentration', fontsize=16)
pyplot.legend(['$N$', '$B_u$', '$B_b$'], fontsize=16)
pyplot.savefig('Freter_wash.png')


# =============================================================================
# fast
# =============================================================================
params=[1, .1, 1, .5, 5, .1, .1, .1, 10, 1, .2, .5]
tspan = np.linspace(0, 4, 500)
N_init=0
Bu_init=.1
Bb_init=0
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [N_init, Bu_init, Bb_init], method='RK45',t_eval=tspan)
pyplot.figure(1)
pyplot.plot(tspan,yp.y[0])
pyplot.plot(tspan,yp.y[1])
pyplot.plot(tspan,yp.y[2])

pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Concentration', fontsize=16)
pyplot.legend(['$N$', '$B_u$', '$B_b$'], fontsize=16)
pyplot.savefig('Freter_fast.png')





# =============================================================================
# Spearman correlation coefficient
# We will also use PyDoe for normal distribution
# =============================================================================



# =============================================================================
#This generates the LHS assuming normal distribution
#==============================================================================
import pyDOE as doe
from scipy.stats.distributions import norm
Num_samples=500
total_bacteria=np.zeros(Num_samples)
parameters=doe.lhs(12,samples=Num_samples)

for i in np.arange(12):
     #Parameters come from normal distribution
     # the 'ppf part provides the values. The standard deviation is the second
     #argument
     parameters[:,i] = norm(loc=params[i], scale=.1*params[i]).ppf(parameters[:, i]) #scaled parameters

N_init=0
Bu_init=.1
Bb_init=0

for i in np.arange(Num_samples):
    params_i=parameters[i,:]
    tspan = np.linspace(0, 100, 500)
    yp= solve_ivp(lambda t,Y: rhs(t,Y,params_i), [tspan[0],tspan[-1]], [N_init, Bu_init, Bb_init ], method='RK45',t_eval=tspan)
    total_bacteria[i]=np.trapz(yp.y[1],tspan)/np.trapz(yp.y[2],tspan)          #trapz is a relatively standard  implementation of the trapezoidal rule for integration
cc=np.zeros([2,12])
cc1=np.zeros([2,12])

for j in np.arange(12):
    cc[:,j]=scipy.stats.pearsonr(parameters[:,j], total_bacteria) 
    
for j in np.arange(12):
    cc1[:,j]=scipy.stats.spearmanr(parameters[:,j], total_bacteria) 



pyplot.figure()
bar_width=.35
pyplot.bar(np.arange(len(Params_names))-.5*bar_width,cc[0], bar_width)
pyplot.bar(np.arange(len(Params_names))+.5*bar_width,cc1[0], bar_width)
pyplot.xticks(np.arange(len(Params_names)), Params_names)
pyplot.legend(['Pearson', 'Spearman'])

#pyplot.savefig('cc_freter.png')

