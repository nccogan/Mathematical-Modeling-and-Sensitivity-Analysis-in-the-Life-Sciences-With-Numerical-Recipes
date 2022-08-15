#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 11:34:10 2021
Rock-Paper-scissors model. DOE model SA based on the paper:
   Sarah C Cotter. A screening design for factorial experiments with
   interactions. Biometrika, 66(2):317–320, 1979. 
@author: cogan
"""



from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.integrate import solve_ivp
import itertools
pyplot.close('all')
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible

#Define the RHS of ODEs
def f1(t,Y,params):
    Pr,Pp,Ps=Y
    alpha_1, alpha_2, alpha_3, beta_1, beta_2, beta_3= params
    return  Pr*(1-Pr-alpha_1*Pp-beta_1*Ps)
def f2(t,Y,params):
    Pr,Pp,Ps=Y
    alpha_1, alpha_2, alpha_3, beta_1, beta_2, beta_3= params
    return  Pp*(1-Pp-alpha_2*Ps-beta_2*Pr)
def f3(t,Y,params):
    Pr,Pp,Ps=Y
    alpha_1, alpha_2, alpha_3, beta_1, beta_2, beta_3= params
    return  Ps*(1-Ps-alpha_1*Pr-beta_1*Pp)

def rhs(t,Y,params):
    Pr,Pp,Ps=Y
    alpha_1, alpha_2, alpha_3, beta_1, beta_2, beta_3= params
    return [f1(t,Y,params),f2(t,Y,params),f3(t,Y,params)]
params_list=[r'$\alpha_1$',r'$\alpha_2$',r'$\alpha_3$', r'$\beta_1$',  r'$\beta_2$',  r'$\beta_3$']
#Change IC to show excitable
params=[.2,.2,.2,2,2,2]
Rr0=.33
Rp0=.1
Rs0=.1
tspan = np.linspace(0, 200, 4000)
y_soln= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [Rr0,Rp0,Rs0], method='LSODA',t_eval=tspan)
pyplot.plot(y_soln.t,y_soln.y[0])
pyplot.plot(y_soln.t,y_soln.y[1])
pyplot.plot(y_soln.t,y_soln.y[2])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Population', fontsize=16)
pyplot.legend(['$P_r$', '$P_p$','$P_s$'], fontsize=16,loc='upper left')
pyplot.savefig('rps_1.png')

pyplot.figure()

params=[.05,.05,.05,3,3,3]
Rr0=.2
Rp0=.1
Rs0=.3
tspan = np.linspace(0, 100, 4000)
y_soln= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [Rr0,Rp0,Rs0], method='LSODA',t_eval=tspan)
pyplot.plot(y_soln.t,y_soln.y[0])
pyplot.plot(y_soln.t,y_soln.y[1])
pyplot.plot(y_soln.t,y_soln.y[2])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Population', fontsize=16)
pyplot.legend(['$P_r$', '$P_p$','$P_s$'], fontsize=16,loc='upper left')
pyplot.savefig('rps_2.png')

# =============================================================================
# make the design matrix indicating high and low values
# =============================================================================


X1=[1,-1] # two factor
number=0
DOE=np.empty([2**len(params),len(params)])
for combination in itertools.product(X1,repeat=len(params)):
    DOE[number,:]=combination
    number+=1
#NOTE that DOE is in backwards order. The last column correspnds to the first parameter etc.
params=[.2,.2,.2,2,2,2]
Rr0=.1
Rp0=.1
Rs0=.1

responser=np.empty([2**len(params)])      #do this for each species
responsep=np.empty([2**len(params)])    
responses=np.empty([2**len(params)])    

for k in np.arange(2**len(params)):
    params1=params+np.multiply(params,.9)*DOE[k,:]
    y_soln= solve_ivp(lambda t,Y: rhs(t,Y,params1), [tspan[0],tspan[-1]], [Rr0,Rp0,Rs0], method='RK45',t_eval=tspan)
    responser[k]=sum(tspan[np.argwhere(y_soln.y[0]>.8)])/tspan[-1]
    responsep[k]=sum(tspan[np.argwhere(y_soln.y[1]>.8)])/tspan[-1]
    responses[k]=sum(tspan[np.argwhere(y_soln.y[2]>.8)])/tspan[-1]

Cotter_SAr=np.empty(len(params))
Cotter_SAp=np.empty(len(params))
Cotter_SAs=np.empty(len(params))

for i in np.arange(len(params)):
    Cotter_SAr[len(params)-1-i]=np.sum(DOE[:,i]*responser)/2**len(params)/2
    Cotter_SAp[len(params)-1-i]=np.sum(DOE[:,i]*responsep)/2**len(params)/2
    Cotter_SAs[len(params)-1-i]=np.sum(DOE[:,i]*responses)/2**len(params)/2
    
bar_width=.3
pyplot.figure()
pyplot.bar(np.arange(len(params_list))-1*bar_width,Cotter_SAr, bar_width)
pyplot.bar(np.arange(len(params_list))-0*.5*bar_width,Cotter_SAp, bar_width)
pyplot.bar(np.arange(len(params_list))+1*bar_width,Cotter_SAs, bar_width)
pyplot.legend(['$P_r$', '$P_p$', '$P_s$'],fontsize=16)
pyplot.xticks(np.arange(len(params_list)), params_list)
pyplot.ylabel('SA',fontsize=16)
pyplot.xlabel('Parameter',fontsize=16)
pyplot.savefig('cotter.png')
