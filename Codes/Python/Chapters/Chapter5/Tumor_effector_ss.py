#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 13:13:54 2021
Slightly different steady-state solver using the scipy package
fsolve.

We also use a generic numerical estimation for the jacobian
based on finite differences. There are more sophisticated ways, but this 
is reasonably effective. 
@author: cogan
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')
#####################from scipy.integrate import ode
import numpy as np
import numpy.matlib
from matplotlib import pyplot    # import plotting library
from scipy.optimize import fsolve
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
params=[ .1181,.3743,1.131,20.19,3.11*10**(-3),1.636,.5*10**3]
#T=Tumor cells
#E=effector cells
#Try and calculate the steady-state. We just need to know reasonable 
#initial guesses.
#We can get a guess for where to start and then use a function to solve
# I am separating the functions out on purpose, I could have used rhs(t,Y,params) to do this as well
#fsolve has a different syntax which
#illustrates one difficulty of using packages -- different input/ouptu structures
##########################
def fs(Y,params):
    (y1,y2) = Y
    f1=params[0]-params[1]*y1+params[2]*y1*y2/(params[3]+y2)-params[4]*y1*y2
    f2=params[5]*y2*(1-y2/params[6])-y1*y2
    return [f1,f2]
#A very simple way to build the Jacobian
def Myjac(t,Y,params,epsilon):
    E,T =Y
    
    f_E = (f(0,[E+epsilon,T],params)-f(0,[E-epsilon,T],params))/epsilon/2
    f_T = (f(0,[E,T+epsilon],params)-f(0,[E,T-epsilon],params))/epsilon/2
    g_E = (g(0,[E+epsilon,T],params)-g(0,[E-epsilon,T],params))/epsilon/2
    g_T = (g(0,[E,T+epsilon],params)-g(0,[E,T-epsilon],params))/epsilon/2
    J = np.array([[f_E,f_T],[g_E,g_T]])
    return J

#Find the steady-states
solution1 = fsolve(fs, (.3,.01),args=(params,) )
solution2 = fsolve(fs, (1,8),args=(params,) )
solution3=fsolve(fs,(.7,280),args=(params,))
solution4=fsolve(fs,(.1,580),args=(params,))



J=Myjac(0,[solution4[0],solution4[1]],params,.001)
w,v= np.linalg.eig(J)
print(w)
