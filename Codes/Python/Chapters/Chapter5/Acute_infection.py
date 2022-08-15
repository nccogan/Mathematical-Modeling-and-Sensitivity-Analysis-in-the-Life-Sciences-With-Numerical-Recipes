#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 22:28:54 2018

@author: cogan

Acute infection model
y1=T: Target cells
y2=I: Infected cells
y3=V: Virus cells
y4=A: antibodies (signals)
params=[beta, delta, p, c, c_a, k, S_a, d]
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
from matplotlib import pyplot    # import plotting library
from scipy.integrate import solve_ivp
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible

pyplot.close('all')

# =============================================================================
# Define the RHS of ODEs
# This section considers the direct
# simulation. 
# 
# =============================================================================
def f1(t,Y,params):
    y1,y2,y3,y4=Y
    return -params[0]*y1*y3
def f2(t,Y,params):
    y1,y2,y3,y4=Y
    return params[0]*y1*y3-params[1]*y2
def f3(t,Y,params):
    y1,y2,y3,y4=Y
    return params[2]*y2-params[3]*y3-params[4]*y3*y4
def f4(t,Y,params):
    y1,y2,y3,y4=Y
    return params[5]*y3*y4+params[6]-params[7]*y4
def rhs(t,Y,params):
    y1,y2,y3,y4=Y
    return [f1(t,Y,params),f2(t,Y,params),f3(t,Y,params),f4(t,Y,params)]


# =============================================================================
# Differential SA
# For this we have augmented the 4 odes with the derivatives with respect to the
#  parameters. We have to indicate whch parameter -- here this is j
# Note that we could be smarter about generalizing our Jacobians rather than typing them by hand
# 
# =============================================================================


def rhs_diff(t,Y,params,j): # j will denote the parameter we are looking at
    y1,y2,y3,y4,s1,s2,s3,s4=Y
    dp=.1*params[j]#*params[j] #increment the parameter
    dy=.001 #for component of Jac
    params_p= np.array(params)
    params_m=np.array(params)
    params_p[j]=params[j]+dp
    params_m[j]=params[j]-dp
#f_p entries
    fp_1 = (f1(t,[y1,y2,y3,y4],params_p)-f1(t,[y1,y2,y3,y4],params_m))/(2*dp)
    fp_2 = (f2(t,[y1,y2,y3,y4],params_p)-f2(t,[y1,y2,y3,y4],params_m))/(2*dp)
    fp_3 = (f3(t,[y1,y2,y3,y4],params_p)-f3(t,[y1,y2,y3,y4],params_m))/(2*dp)
    fp_4 = (f4(t,[y1,y2,y3,y4],params_p)-f4(t,[y1,y2,y3,y4],params_m))/(2*dp)
#We could loop through to get the Jacobian but we are doing it explicitly
#First row
    f1_1 = (f1(t,[y1+dy,y2,y3,y4],params)-f1(t,[y1-dy,y2,y3,y4],params))/(2*dy)
    f1_2 = (f1(t,[y1,y2+dy,y3,y4],params)-f1(t,[y1,y2-dy,y3,y4],params))/(2*dy)
    f1_3 = (f1(t,[y1,y2,y3+dy,y4],params)-f1(t,[y1,y2,y3-dy,y4],params))/(2*dy)
    f1_4 = (f1(t,[y1,y2,y3,y4+dy],params)-f1(t,[y1,y2,y3,y4-dy],params))/(2*dy)
#Second row
    f2_1 = (f2(t,[y1+dy,y2,y3,y4],params)-f2(t,[y1-dy,y2,y3,y4],params))/(2*dy)
    f2_2 = (f2(t,[y1,y2+dy,y3,y4],params)-f2(t,[y1,y2-dy,y3,y4],params))/(2*dy)
    f2_3 = (f2(t,[y1,y2,y3+dy,y4],params)-f2(t,[y1,y2,y3-dy,y4],params))/(2*dy)
    f2_4 = (f2(t,[y1,y2,y3,y4+dy],params)-f2(t,[y1,y2,y3,y4-dy],params))/(2*dy)
#Third row
    f3_1 = (f3(t,[y1+dy,y2,y3,y4],params)-f3(t,[y1-dy,y2,y3,y4],params))/(2*dy)
    f3_2 = (f3(t,[y1,y2+dy,y3,y4],params)-f3(t,[y1,y2-dy,y3,y4],params))/(2*dy)
    f3_3 = (f3(t,[y1,y2,y3+dy,y4],params)-f3(t,[y1,y2,y3-dy,y4],params))/(2*dy)
    f3_4 = (f3(t,[y1,y2,y3,y4+dy],params)-f3(t,[y1,y2,y3,y4-dy],params))/(2*dy)
#Second row
    f4_1 = (f4(t,[y1+dy,y2,y3,y4],params)-f4(t,[y1-dy,y2,y3,y4],params))/(2*dy)
    f4_2 = (f4(t,[y1,y2+dy,y3,y4],params)-f4(t,[y1,y2-dy,y3,y4],params))/(2*dy)
    f4_3 = (f4(t,[y1,y2,y3+dy,y4],params)-f4(t,[y1,y2,y3-dy,y4],params))/(2*dy)
    f4_4 = (f4(t,[y1,y2,y3,y4+dy],params)-f4(t,[y1,y2,y3,y4-dy],params))/(2*dy)
    return [f1(t,[y1,y2,y3,y4],params),f2(t,[y1,y2,y3,y4],params),f3(t,[y1,y2,y3,y4],params),f4(t,[y1,y2,y3,y4],params),
            fp_1+f1_1*s1+f1_2*s2+f1_3*s3+f1_4*s4,
            fp_2+f2_1*s1+f2_2*s2+f2_3*s3+f2_4*s4,
            fp_3+f3_1*s1+f3_2*s2+f3_3*s3+f3_4*s4,
            fp_4+f4_1*s1+f4_2*s2+f4_3*s3+f4_4*s4]            

# =============================================================================
# Direct numerical simulation
# =============================================================================

params=[ .1181,.0743,1.131,20.19,3.11*10**(-3),1.636,.5*10**3,1]

tspan = np.linspace(0, 150, 500)

y0 = [100.35,1,0,0]
y_simulation = solve_ivp(lambda t,Y: rhs(t,Y,params),  [tspan[0],tspan[-1]],y0,method='RK45',t_eval=tspan)
 
pyplot.figure(1)
pyplot.plot(tspan,y_simulation.y[0])
pyplot.plot(tspan,y_simulation.y[1])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Concentration', fontsize=16)
pyplot.legend(['Target','Infected'], fontsize=16)

pyplot.figure(2)
pyplot.plot(tspan,y_simulation.y[2])
pyplot.legend(['Virus'], fontsize=16)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Viral Concentration', fontsize=16)


pyplot.figure(3)
pyplot.plot(tspan,y_simulation.y[3])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Antibody Concentration', fontsize=16)


#Post-processing to estimate derivative using the gradient operator
pyplot.figure(4)
gradient_v=np.gradient(y_simulation.y[3],(tspan[2]-tspan[1]));
pyplot.plot(tspan, gradient_v)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel(r'$\frac{dV}{dt}$', fontsize=16)


# =============================================================================
# Feature Sensitivity using post-processing
# =============================================================================
tspan = np.linspace(0, 100, 100)
y_SA_feature = solve_ivp(lambda t,Y: rhs_diff(t,Y,params,6),  [tspan[0],tspan[-1]],[100.35,1,0,0,0,0,0,0],method='RK45',t_eval=tspan)

pyplot.close('all')

pyplot.figure(1)
pyplot.plot(tspan,y_SA_feature.y[2])
pyplot.figure(2)
pyplot.plot(tspan,y_SA_feature.y[6])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Sens. of $V$ to $S_A$', fontsize=16)


gradient_v_feature=np.gradient(y_SA_feature.y[6],(tspan[2]-tspan[1]));

pyplot.figure(4)
pyplot.plot(tspan, gradient_v_feature)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Sens. of ' r'$\frac{\partial V}{\partial t}$ to $S_A$', fontsize=16)





