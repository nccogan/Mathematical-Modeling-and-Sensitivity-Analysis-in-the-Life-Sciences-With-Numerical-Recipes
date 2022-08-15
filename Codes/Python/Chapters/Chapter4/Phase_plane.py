"""
Created on Thu Sep 20 12:02:46 2018

@author: cogan
Adapted from 
http://kitchingroup.cheme.cmu.edu/blog/2013/02/21/Phase-portraits-of-a-system-of-ODEs/

"""
from IPython import get_ipython
get_ipython().magic('reset -sf')

import numpy as np
import matplotlib.pyplot as pyplot
from scipy.integrate import solve_ivp

pyplot.close('all')

def f(t,Y,params):
    y1,y2=Y
    alpha, beta, delta, gamma=params
    return alpha*y1*y2-delta*y1
def g(t,Y,params):
    y1,y2=Y
    alpha, beta, delta, gamma=params
    return -beta*y1*y2+gamma*y2*(10-y2)
def rhs(t,Y,params):
    y1,y2=Y
    alpha, beta, delta, gamma=params
    return [f(t,Y,params),g(t,Y,params)]

#Define the parameters    
#params=[1,1,.2, 1, 1, 1] 
params=[.2,.1,.2,.1]
#Define the time to simulate
Tfinal=200
#Define the time discretization
tspan = np.linspace(0, Tfinal, 500)
#Define Initial Conditions
IC = [.1, 1]
#Solve the ODE
yp = solve_ivp(lambda t,Y: rhs(t,Y,params), [0,Tfinal], IC, method='RK45',t_eval=tspan)

#Phase-plane View
      
y1 = np.arange(0, 22.0,1)
y2 = np.arange(0, 10.0,.5)
Y1, Y2 = np.meshgrid(y1, y2)

t = 0

pyplot.quiver(Y1, Y2, f(0,[Y1,Y2],params),g(0,[Y1,Y2],params),color='#FF800E')
pyplot.plot(yp.y[0], yp.y[1]) # start
pyplot.xlabel('$y_1$')
pyplot.ylabel('$y_2$')
