
"""
Created on Wed Dec 22 21:11:32 2021 Basic ODE Script
@author: cogan
"""
from IPython import get_ipython 
get_ipython().magic('reset -sf')
import numpy as np
from matplotlib import pyplot
# import plotting library
from scipy.integrate import solve_ivp
#Define the right-hand-side of the differential equation
pyplot.close('all')
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible

def rhs(t,Y,params): 
    r, K = params 
    y=Y
    return [r*y*(K-y)]
 # Initial Condition
IC=.2
#Intial time
tstart=0
#Final time
tstop=4
tspan = np.linspace(tstart, tstop, 100) 
#Parameters
params=[.1,100]
y_solution = solve_ivp(lambda t,Y: rhs(t,Y,params), [tstart,tstop], [IC], method='RK45',t_eval=tspan)
pyplot.plot(y_solution.t, y_solution.y[0],'k',label='y(t)') 
pyplot.ylabel('y(t) ',fontsize=16) 
pyplot.xlabel('t',fontsize=16) 
#pyplot.savefig('Logistic.png')
