#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 11:11:10 2022
Differential Sensitivity for the 
Tumor/Effector system
@author: cogan
"""

from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import numpy.matlib
from matplotlib import pyplot    # import plotting library
from scipy.integrate import solve_ivp
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible

pyplot.close('all')


# same as usual: define f and g, rhs returns the right hand side
# parameters are from the paper, beware the order though

def f(t,Y,params):
    E,T =Y
    s,d,p,g,m,r,k=params
    return s-d*E+p*E*T/(g+T)-m*E*T
def g(t,Y,params):
    E,T =Y
    s,d,p,g,m,r,k=params
    return r*T*(1-T/k)-E*T
def rhs(t,Y,params):
    E,T =Y
    s,d,p,g,m,r,k=params
    return [f(t,Y,params),g(t,Y,params)]

def rhs_diff(t,Y,params,j): # j will denote the parameter we are looking at
    e,t,s1,s2=Y
    dp=.1*params[j]#*params[j] #increment the parameter
    de=.001 #for e component of Jac
    dt=.001 #for t component of Jac
    params_p= np.array(params)
    params_m=np.array(params)
    params_p[j]=params[j]+dp
    params_m[j]=params[j]-dp
    f1 = (f(0,[e,t],params_p)-f(0,[e,t],params_m))/(2*dp)
    g1 = (g(0,[e,t],params_p)-g(0,[e,t],params_m))/(2*dp)
    f_e = (f(0,[e+de,t],params)-f(0,[e-de,t],params))/(2*de)
    f_t = (f(0,[e,t+dt],params)-f(0,[e,t-dt],params))/(2*dt)
    g_e = (g(0,[e+de,t],params)-g(0,[e-de,t],params))/(2*de)
    g_t = (g(0,[e,t+dt],params)-g(0,[e,t-dt],params))/(2*dt)
    #[f_e,f_t;g_e,g_t] is the Jac
    return [f(0,[e,t],params),g(0,[e,t],params),f1+f_e*s1+f_t*s2,g1+g_e*s1+g_t*s2] 
params=[ .1181,.3743,1.131,20.19,3.11*10**(-3),1.636,.5*10**3]
#T=Tumor cells
#E=effector cells 
tspan = np.linspace(0, 100.000, 200)
# =============================================================================
# There are many different sensitivies we can explore
#There are four different steady-states and 7 parameters
#Here SA_1, SA_2, SA_3, and SA_4 denote the sensitivities near
#different steady-states. Each for loop calculates the solution
#to the augmented system
# =============================================================================
SA_1=np.empty((7,2))
for j in np.arange(0,7):
    y_SA1 = solve_ivp(lambda t,Y: rhs_diff(t,Y,params,j),  [tspan[0],tspan[-1]],[.1,.10,0,0],method='LSODA',t_eval=tspan)
    SA_1[j,:]=y_SA1.y[-1,2:4]

SA_2=np.empty((7,2))
for j in np.arange(0,7):
    y_SA2 = solve_ivp(lambda t,Y: rhs_diff(t,Y,params,j),  [tspan[0],tspan[-1]],[.1,10,0,0],method='LSODA',t_eval=tspan)
    SA_2[j,:]=y_SA2.y[-1,2:4]

SA_3=np.empty((7,2))
for j in np.arange(0,7):
    y_SA3 = solve_ivp(lambda t,Y: rhs_diff(t,Y,params,j),  [tspan[0],tspan[-1]],[1, 280,0,0],method='LSODA',t_eval=tspan)
    SA_3[j,:]=y_SA3.y[-1,2:4]

SA_4=np.empty((7,2))
for j in np.arange(0,7):
    y_SA4 = solve_ivp(lambda t,Y: rhs_diff(t,Y,params,j),  [tspan[0],tspan[-1]],[0.001, 440 ,0,0],method='LSODA',t_eval=tspan)
    SA_4[j,:]=y_SA4.y[-1,2:4]
    
# =============================================================================
# These can be plotted using a bar plot to compare
#how the tumor and effector QoIs compare
# =============================================================================
p_names=('s', 'd', 'p', 'g', 'm', 'r', 'k')
bar_width = 0.35
pyplot.bar(np.arange(len(p_names))-bar_width/2,SA_4[:,0], bar_width)
pyplot.bar(np.arange(len(p_names))+bar_width/2,SA_4[:,1], bar_width)
pyplot.legend(['$S_{E,i}$', '$S_{T,i}$'],fontsize=16)
pyplot.xticks(np.arange(len(p_names)), p_names)

# =============================================================================
# To investigate how the direct sensitivities change in time for each parameter
#we can look at the solution to the system
# =============================================================================
pyplot.figure()
j=0 #which parameter
y_time = solve_ivp(lambda t,Y: rhs_diff(t,Y,params,j),  [tspan[0],tspan[-1]],[1, 10 ,0,0],method='LSODA',t_eval=tspan)
pyplot.plot(y_time.t, y_time.y[2]) #The sensitivity of the Effector to the jth paramter
pyplot.plot(y_time.t, y_time.y[3]) #The sensitivity of the Tumor to the jth paramter
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Sensitivities', fontsize=16)
pyplot.legend(['Tumor','Effector'])

