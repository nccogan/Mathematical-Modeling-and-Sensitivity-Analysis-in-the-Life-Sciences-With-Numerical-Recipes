#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:31:09 2022

@author: cogan
"""

"""
Created on Thu Dec 30 17:13:57 2021
Competitive Exclusion
With random sampling of parameters and linear regression 
to approximate the slope
Note that the closer the relationship is to linear, 
the more appropriate the conclusions are. 
@author: cogan
"""

from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
from matplotlib import pyplot    # import plotting library
from scipy.integrate import solve_ivp
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible
pyplot.close('all')

def f(t,Y,params):
    y1,y2=Y
    r1, kappa1, alpha12, r2, kappa2, alpha21 = params
    return r1*y1*(kappa1-y1-alpha12*y2)/kappa1
def g(t,Y,params):
    y1,y2=Y
    r1, kappa1, alpha12, r2, kappa2, alpha21 = params    
    return r2*y2*(kappa2-y2-alpha21*y1)/kappa2
def rhs(t,Y,params):
    y1,y2=Y
    r1, kappa1, alpha12, r2, kappa2, alpha21 = params
    return [f(t,Y,params),g(t,Y,params)]

#Parameter Descriptions
#params=[r1,kappa1,alpha12,r2,kappa2,alpha21]
#Example parameters for different cases of the
#competitive exclusion
#case 1 params=[1,1,2,1,2,1]
#case 2 
params=[1,2,1,1,1,2]
#case 3 params=[1,2,1,1,3,2]
#case 4 params=[1,3,2,1,2,1]



#Sample parameters
Num_samples=500
tspan=np.linspace(0,26,200)
#Initial conditions for N1 and N2
N1_0=.1
N2_0=.1
y_solution = solve_ivp(lambda t,Y: rhs(t,Y,params),  [tspan[0],tspan[-1]],[N1_0, N2_0],method='RK45',t_eval=tspan)

pyplot.plot(tspan,y_solution.y[0])
pyplot.plot(tspan,y_solution.y[1])
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Population', fontsize=16)

pyplot.legend(['$N_1$','$N_2$'], fontsize=16)
#pyplot.savefig('CE_State4.png')
