#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 12:49:41 2018
Fixed point (Steady-State) solver, example
@author: cogan
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')
#####################from scipy.integrate import ode
import numpy as np
import numpy.matlib
from matplotlib import pyplot    # import plotting library


fp = []
def f(t,Y,params):
    x,y=Y
    alpha, beta, delta, gamma =params
    return alpha*x*y- delta* x

def g(t,Y,params):
    x,y=Y
    alpha, beta, delta, gamma =params
    return -beta*x*y+gamma*y

#loops throught the given ranges of x and y to evaluate whether
#f(x,y)=0 and g(x,y)=0
def find_fixed_points(r):
    for x in range(r):
        for y in range(r):
            t=[0]
            if ((f(t,[x,y],params) == 0) and (g(t,[x,y],params) == 0)):
                fp.append((x,y))
                print('The system has a fixed point  %s,%s' % (x,y))
    return fp

params=[.2,.1,.2,.1]
find_fixed_points(2)


#Contour plot
x = np.arange(-2,2,.1)
y = np.arange(-2,2,.1)
X, Y = np.meshgrid(x, y)
pyplot.contour(X,Y,f(0,[X,Y],params),colors='C1',levels=0)
pyplot.contour(X,Y,g(0,[X,Y],params),colors='C2',levels=0)
