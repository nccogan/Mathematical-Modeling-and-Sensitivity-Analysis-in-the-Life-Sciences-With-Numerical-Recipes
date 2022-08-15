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
#case 2 params=[1,2,1,1,1,2]
#case 3
params=[1,2,1,1,3,2]
#case 4params=[1,3,2,1,2,1]



#Sample parameters
Num_samples=500
tstart=0
tstop=26
#Initial conditions for N1 and N2
N1_0=.1
N2_0=.1
#Bound the parameter sample space
params_max=np.multiply(1.4,[1,2,1,1,3,2])
params_min=np.multiply(.6,[1,2,1,1,3,2])
#just for labelling
pv=['$r_1$','$\\kappa_1$','$\\alpha_{12}$','$r_2$','$\\kappa_2$','$\\alpha_{21}$']
#Initialize QoI
p_vals=np.empty(shape=Num_samples)
#Initialize different QoIs
QoI1=np.empty(shape=(Num_samples))
QoI2=np.empty(shape=(Num_samples))


#This loop runs through all parameters. 
#The outside loop (k)
#fixes the parameter of interest.
#The inside loop (s) 
#goes through each sample

for k in np.arange(0,6): 
    params=(params_max+params_min)/2
    for s in np.arange(0,Num_samples):
        params[k]=np.random.uniform(params_min[k],params_max[k], size=None)
        tspan = np.linspace(tstart, tstop, 200)
        Y0 = [N1_0, N2_0]
        y_solution = solve_ivp(lambda t,Y: rhs(t,Y,params),  [tspan[0],tspan[-1]],Y0,method='RK45',t_eval=tspan)
        ys=y_solution.y
 #       solve_ivp(lambda t,Y: f(t,Y,params),  [tspan[0],tspan[-1]],y0,method='RK45',t_eval=tspan)
        QoI1[s]=ys[-1,0]/(ys[-1,0]+ys[-1,1])
        QoI2[s]=ys[-1,1]/(ys[-1,0]+ys[-1,1])
        p_vals[s]=params[k]
    pyplot.scatter(p_vals, QoI1)
    pyplot.ylabel('$QoI1$', fontsize = 16)
    pyplot.xlabel('%s'%pv[k], fontsize = 16)
    coeffs=np.polyfit(p_vals,QoI1,  1)  #Polyfit interpolates the particular QoI vs parameter data 
                                        #generated with a polynomial. The degree of the polynomial is an option an the return is the 
                                        #coefficent. 
    pyplot.figtext(.3,.2,'Slope  = %10f' %coeffs[0], fontsize = 16)
    pyplot.plot(p_vals,coeffs[1]+coeffs[0]*p_vals,color='#C85200',linewidth=3)
    pyplot.figure()
#    pyplot.savefig('reg_final%d.png' %k,bbox_inches="tight")
 #   pyplot.close('all')
