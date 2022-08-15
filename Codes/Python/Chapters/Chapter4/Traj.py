#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 10:46:44 2018

@author: cogan
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 12:02:46 2018

@author: cogan

"""

from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.integrate import solve_ivp
pyplot.close('all')
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible


#Define the RHS, FIrst for Lotka-Volterra
def f(t,Y,params):
    y=Y
    r,k=params
    return r*y*(1-y/k)
        
params=[.2,10]

t = 0

tspan = np.linspace(0, 50, 500)
y0 = [.1]
ys= solve_ivp(lambda t,Y: f(t,Y,params), [tspan[0],tspan[-1]], y0, method='RK45',t_eval=tspan)

pyplot.plot(tspan, ys.y[0]) # start

pyplot.xlim([tspan[0],tspan[-1]])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Population', fontsize=16)

#pyplot.savefig('logistic.png')
