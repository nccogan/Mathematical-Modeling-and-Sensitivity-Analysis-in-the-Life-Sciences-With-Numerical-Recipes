#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 31 22:28:06 2021
TB Model
@author: cogan
B= Free bacteria
Q = Dormant Bacteria
X= Immune response
r,g,h,f,f,g,a,s,k,d=params
"""

from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
from matplotlib import pyplot    # import plotting library
from scipy.integrate import solve_ivp
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible


#Define the RHS, FIrst for Lotka-Volterra
def B_rhs(t,Y,params):
    B,Q,X=Y
    r,h,f,g,a,s,k,d=params
    return r*B+g*Q-(h*B*X+f*B)
def Q_rhs(t,Y,params):
    B,Q,X=Y
    r,h,f,g,a,s,k,d=params
    return f*B-g*Q
def X_rhs(t,Y,params):
    B,Q,X=Y
    r,h,f,g,a,s,k,d=params
    return a+s*X*(B/(k+B))-d*X
                  
def rhs(t,Y,params):
    B,Q,X=Y
    r,h,f,g,a,s,k,d=params
    return [B_rhs(t,Y,params),Q_rhs(t,Y,params),X_rhs(t,Y,params)]

params=[1,.001,.5,.1,.1,1,100,.1]

tspan = np.linspace(0, 160, 200)
y0 = [1,1,1]
y_solution = solve_ivp(lambda t,Y: rhs(t,Y,params),  [tspan[0],tspan[-1]],y0,method='RK45',t_eval=tspan)
y_plot=y_solution.y  #take the y values

pyplot.figure()
pyplot.plot(tspan, np.log(y_plot[0]))
pyplot.plot(tspan, np.log(y_plot[0]+y_plot[1]))
pyplot.plot(tspan, np.log(y_plot[2]))
pyplot.legend(['B','B+Q','X'],fontsize=16)
pyplot.xlabel('Time',fontsize=16)
pyplot.ylabel('log(density)',fontsize=16)
pyplot.savefig('tb.png')

##########################################
#function to calculate the timeof maximum
######################################


Num_samples=100
params_max=np.multiply(params,1.1)
params_min=np.multiply(params,.9)
time_max=np.empty(shape=Num_samples)
for s in np.arange(0,Num_samples):   
    for k in np.arange(0,6): 
        params[k]=np.random.uniform(params_min[k],params_max[k], size=None)
    y_solution = solve_ivp(lambda t,Y: rhs(t,Y,params),  [tspan[0],tspan[-1]],y0,method='RK45',t_eval=tspan)
    time_max[s]=tspan[np.argmax(y_solution.y[0])]
pyplot.figure()

pyplot.hist(time_max,15)
pyplot.xlabel('Time of Maximum',fontsize=16)
pyplot.ylabel('Frequency',fontsize=16)
pyplot.savefig('time_max_hist.png')

    

##################
#Local Sensitivity
######################################
Num_samples=10
params_max=np.multiply(params,1.2)
params_min=np.multiply(params,.9)
time_max=np.empty(shape=Num_samples)
s=np.empty(shape=7) #place holder for der. of QoI
r=np.empty(shape=7) #place holder for params

time_max_min=np.empty(shape=7)
time_max_max=np.empty(shape=7)
params_small=[1,.001,.5,.1,.1,1,100,.1]
params_large=[1,.001,.5,.1,.1,1,100,.1]

for k in np.arange(0,7): 
    params_small[k]=params_min[k]
    y_solution_min = solve_ivp(lambda t,Y: rhs(t,Y,params_small),  [tspan[0],tspan[-1]],y0,method='RK45',t_eval=tspan)
    time_max_min[k]=tspan[np.argmax(y_solution_min.y[0])]
    params_large[k]=params_max[k]
    y_solution_max = solve_ivp(lambda t,Y: rhs(t,Y,params_large),  [tspan[0],tspan[-1]],y0,method='RK45',t_eval=tspan)
    time_max_max[k]=tspan[np.argmax(y_solution_max.y[0])]
    s[k]=(time_max_max[k]-time_max_min[k])/tspan[np.argmax(y_solution_max.y[0])]
    r[k]=(params_large[k]-params_small[k])/params[k]
    params=[1,.001,.5,.1,.1,1,100,.1]
    params_small=[1,.001,.5,.1,.1,1,100,.1]
    params_large=[1,.001,.5,.1,.1,1,100,.1]
pyplot.figure()
bar_width=.35
p_names=('$r$', '$g$', '$h$', '$f$', '$a$','$s$','$k$' )
pyplot.bar(np.arange(len(p_names)),s/r, bar_width)
pyplot.xticks(np.arange(len(p_names)), p_names)
pyplot.savefig('sa_bar.png')




















