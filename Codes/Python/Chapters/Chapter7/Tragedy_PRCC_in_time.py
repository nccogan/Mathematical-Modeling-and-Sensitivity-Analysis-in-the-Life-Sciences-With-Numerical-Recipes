
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 11:34:10 2021
Cheater model. PRCC in time is implemented with pyDOE
@author: cogan
"""



from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.integrate import solve_ivp
import pandas as pd
import pyDOE as doe
import pingouin as pg
from scipy.stats.distributions import norm

 #   import pingouin as pg

pyplot.close('all')
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible


#Define the RHS of ODEs
def f1(t,Y,params): #S
    S, P, E, B1, B2=Y
    F, S0, gamma, Yield1, mu1, K1, Yield2, mu2, K2, alpha = params
    return F*S0-F*S-gamma*S*E
def f2(t,Y,params): #P
    S, P, E, B1, B2=Y
    F, S0, gamma, Yield1, mu1, K1, Yield2, mu2, K2, alpha = params
    return gamma*S*E-1/Yield1*mu1*P/(K1+P)*B1-1/Yield2*mu2*P/(K2+P)*B2-F*P
def f3(t,Y,params): #E
    S, P, E, B1, B2=Y
    F, S0, gamma, Yield1, mu1, K1, Yield2, mu2, K2, alpha = params
    return alpha*mu1*P*B1/(K1+P)-F*E
def f4(t,Y,params): # B1
    S, P, E, B1, B2=Y
    F, S0, gamma, Yield1, mu1, K1, Yield2, mu2, K2, alpha = params
    return (1-alpha)*mu1*P*B1/(K1+P)-F*B1
def f5(t,Y,params): #B2
    S, P, E, B1, B2=Y
    F, S0, gamma, Yield1, mu1, K1, Yield2, mu2, K2, alpha = params
    return mu2*P*B2/(K2+P)-F*B2
        
def rhs(t,Y,params):
    S, P, E, B1, B2=Y
    F, S0, gamma, Yield1, mu1, K1, Yield2, mu2, K2, alpha = params
    return [f1(t,Y,params),f2(t,Y,params),f3(t,Y,params),f4(t,Y,params),f5(t,Y,params)]

Params_names=['$F$', '$S_0$', r'$\gamma$', '$Y_1$', r'$\mu_1$', '$K_1$', '$Y_2$', r'$\mu_2$', '$K_2$',  r'$\alpha$'] 
# =============================================================================
# Note the differing initial conditions, here cheaters are present
# =============================================================================
params=[1, 1, 20,1, 5, .05, 1, 5, .05, .2]
tspan = np.linspace(0, 25, 500)
S_init=1
P_init=0
E_init=0.1
B1_init=.2
B2_init=0.02
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [S_init, P_init, E_init, B1_init, B2_init], method='LSODA',t_eval=tspan)
pyplot.figure()

pyplot.plot(tspan,yp.y[0])
pyplot.plot(tspan,yp.y[1])
pyplot.plot(tspan,yp.y[2])
pyplot.plot(tspan,yp.y[3])
pyplot.plot(tspan,yp.y[4])

pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Concentration', fontsize=16)
pyplot.legend(['$S$', '$P$', '$E$', r'$B_1$', '$B_2$'], fontsize=16)
pyplot.title(r'$\mu_1=5$ and $\bar{S}$ = 0.02 ' , fontsize=16)

# =============================================================================
# Note the differing initial conditions, here cheaters are absent
# =============================================================================

params=[1, 1, 20,1, 5, .05, 1, 5, .05, .2]
tspan = np.linspace(0, 25, 500)
S_init=1
P_init=0
E_init=0.1
B1_init=.2
B2_init=0*0.02
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [S_init, P_init, E_init, B1_init, B2_init], method='LSODA',t_eval=tspan)
pyplot.figure()


pyplot.plot(tspan,yp.y[0])
pyplot.plot(tspan,yp.y[1])
pyplot.plot(tspan,yp.y[2])
pyplot.plot(tspan,yp.y[3])
pyplot.plot(tspan,yp.y[4])

pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Concentration', fontsize=16)
pyplot.title(r'$\mu_1=5$ and $\bar{S}$ = 0 ' , fontsize=16)
pyplot.legend(['$S$', '$P$', '$E$', r'$B_1$', '$B_2$'], fontsize=16)

# =============================================================================
# Note the differing parameters, here mu_1 is small, no cheaters but the system 
# washout case is stable
# =============================================================================

params=[1, 1, 20,1, .5, .05, 1, 5, .05, .2]
tspan = np.linspace(0, 25, 500)
S_init=1
P_init=0
E_init=0.1
B1_init=.2
B2_init=0*0.02
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [S_init, P_init, E_init, B1_init, B2_init], method='LSODA',t_eval=tspan)
pyplot.figure()

pyplot.plot(tspan,yp.y[0])
pyplot.plot(tspan,yp.y[1])
pyplot.plot(tspan,yp.y[2])
pyplot.plot(tspan,yp.y[3])
pyplot.plot(tspan,yp.y[4])

pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Concentration', fontsize=16)
pyplot.legend(['$S$', '$P$', '$E$', r'$B_1$', '$B_2$'], fontsize=16)
pyplot.title(r'$\mu_1=.5$ and $\bar{S}$ = 0 ' , fontsize=16)



# =============================================================================
# Note the differing parameters, here mu_1 is small, no cheaters but the system 
# washout case is almost stable
# =============================================================================

params=[1, 1, 20,1, 1.45, .05, 1, 5, .05, .2]
tspan = np.linspace(0, 25, 500)
S_init=1
P_init=0
E_init=0.1
B1_init=.2
B2_init=0*0.02
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [S_init, P_init, E_init, B1_init, B2_init], method='LSODA',t_eval=tspan)
pyplot.figure()

pyplot.plot(tspan,yp.y[0])
pyplot.plot(tspan,yp.y[1])
pyplot.plot(tspan,yp.y[2])
pyplot.plot(tspan,yp.y[3])
pyplot.plot(tspan,yp.y[4])

pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Concentration', fontsize=16)
pyplot.legend(['$S$', '$P$', '$E$', r'$B_1$', '$B_2$'], fontsize=16)
pyplot.title(r'$\mu_1=1.45$ and $\bar{S}$ = 0 ' , fontsize=16)


# =============================================================================
#partial correlation coefficient
#We will also use PyDoe for normal distribution
# =============================================================================

# =============================================================================
#This generates the LHS assuming normal distribution
#==============================================================================
import pyDOE as doe
from scipy.stats.distributions import norm
Num_samples=50
total_bacteria=np.zeros(Num_samples)
parameters=doe.lhs(10,samples=Num_samples)

for i in np.arange(10):
     parameters[:,i] = norm(loc=params[i], scale=.025*params[i]).ppf(parameters[:, i]) #scaled parameters


params=[1, 1, 20,1, 1.45, .05, 1, 5, .05, .2]
tspan = np.linspace(0, 25, 500)
S_init=1
P_init=0
E_init=0.1
B1_init=.2
B2_init=0.02
PRCC_in_time=np.empty([11,10])
# =============================================================================
# Time loop for moving the ending time out time: we are using LSODA as the integrator
# =============================================================================
for k in np.arange(9):
    
    tfinal=[.001,.01,.1,1,5,10,30, 40, 50,100][k]
    correlation_matrix=np.empty([Num_samples,11])
    for i in np.arange(Num_samples):
        params_i=parameters[i,:]
        tspan = np.linspace(0, tfinal, 500)
        yp= solve_ivp(lambda t,Y: rhs(t,Y,params_i), [tspan[0],tspan[-1]], [S_init, P_init, E_init, B1_init, B2_init], method='LSODA',t_eval=tspan)
        correlation_matrix[i,:]=np.append(parameters[i,:],yp.y[3][-1])

    Correlation_data=pd.DataFrame(data=correlation_matrix,
              index=pd.RangeIndex(range(0, Num_samples)),
              columns=pd.RangeIndex(range(0, 11)))
    Correlation_data.columns = ['$F$', '$S_0$', r'$\gamma$', '$Y_1$', r'$\mu_1$', '$K_1$', '$Y_2$', r'$\mu_2$', '$K_2$',  r'$\alpha$', 'QoI']
    
    PRCC=Correlation_data.pcorr()
    PRCC_in_time[:,k]=PRCC[PRCC.columns[-1]].values
 
for k in np.arange(10):
    pyplot.plot([.001,.01,.1,1,5,10,30, 40],PRCC_in_time[k,0:-2],'o-')
    pyplot.legend(['$F$', '$S_0$', r'$\gamma$', '$Y_1$', r'$\mu_1$', '$K_1$', '$Y_2$', r'$\mu_2$', '$K_2$',  r'$\alpha$'], bbox_to_anchor=(1.01,1), loc='upper left', fontsize=16)
#    pyplot.xticks(np.arange(7),['T=.001','T=..01','T=..1','T=.1','T=.5','T=.10','T=.50'])
    pyplot.xlabel(r'$t_{final}$', fontsize=16)
    pyplot.ylabel(r'$PRCC$', fontsize=16)



