#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 13:46:41 2019

@author: cogan
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 10:55:12 2019

@author: cogan
Coding Pearson (non-ranked transform) and Spearman
the main idea is to make a matrix of inputs (parameter sets)
and column of outputs (model evaluated at those parameters) 
and consider the correlation between the columns of parameters and column of outputs.
Recalling that the column of parameters is just different realizations of a single parameter
P = [X_p1;X_p2;X_p3...], X_pi is a column vector of parameters. 
Pearson: Cov(X_pi,Y_out)/sqrt(Var(X_pi)*var(y))
@author: cogan
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import numpy.matlib
from scipy.integrate import odeint
import matplotlib.pyplot as pyplot
from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.integrate import solve_ivp
import scipy as scipy
pyplot.close('all')
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible

#Define the RHS of ODEs
def f1(t,Y,params):
    S,P,E,X1,X2=Y
    S0,D,k1,k2,Y,mu,K,q=params
    return D*(S0-S)-(k1*S*E-k2*P)
def f2(t,Y,params):
    S,P,E,X1,X2=Y
    S0,D,k1,k2,Y,mu,K,q=params
    return (k1*S*E-k2*P)-1/Y*(X1+X2)*mu*P/(K+P)-D*P
def f3(t,Y,params):
    S,P,E,X1,X2=Y
    S0,D,k1,k2,Y,mu,K,q=params
    return (1-q)*X1*mu*P/(K+P)-D*E
def f4(t,Y,params):
    S,P,E,X1,X2=Y
    S0,D,k1,k2,Y,mu,K,q=params
    return X1*q*mu*P/(K+P)-D*X1
def f5(t,Y,params):
    S,P,E,X1,X2=Y
    S0,D,k1,k2,Y,mu,K,q=params
    return X2*mu*P/(K+P)-D*X2
def rhs(t,Y,params):
    S,P,E,X1,X2=Y
#    S0,D,k1,k2,Y,mu,K,q=params
    return [f1(t,Y,params),f2(t,Y,params),f3(t,Y,params),f4(t,Y,params),f5(t,Y,params)]
#####################from scipy.integrate import ode
param_list=['S','D','k1','k2','Y','mu','K','q']
params_0=[ 1,1, 20, .005, 1, 5,.05,.8] #Nominal parameter values
tspan1 = np.linspace(0, 25, 500)
y0=[1,0,.1,.2,.1]
yp= solve_ivp(lambda t,Y: rhs(t,Y,params_0), [tspan1[0],tspan1[-1]], y0, method='RK45',t_eval=tspan1)
pyplot.figure(1) 
pyplot.plot(tspan1,yp.y[0])
pyplot.plot(tspan1,yp.y[1])
pyplot.plot(tspan1,yp.y[2])
pyplot.plot(tspan1,yp.y[3])
pyplot.plot(tspan1,yp.y[4])

# =============================================================================
# Make the parameter matrix, we could 
# =============================================================================
N=1000 # number of samples
del_p=.25 # parameter variation e.g. 10%
p_len=len(params_0)
X_param=np.zeros((p_len,N))
Y_out=np.zeros((N))
SA_pearson=np.zeros(p_len)
SA_spearman=np.zeros(p_len)
SA_PRCC=np.zeros(p_len)

tspan = np.linspace(0, 3, 500)

for i in range(0,p_len):
    #one parameter
    X_param[i,:]=np.random.uniform(params_0[i]*(1-del_p),params_0[i]*(1+del_p),N) 

for k in range(0,N):
    params=X_param[:,k]  
    yp = odeint(rhs, y0,tspan,args=(params,))
    Y_out[k]=yp[3,-1]+yp[4,-1]#max(yp[1,:])
#    Y_out[k]=max(yp[1,:])
v_y=np.var(Y_out, dtype=np.float64)

for k in range(0,N):
    params=X_param[:,k]  
    yp = odeint(rhs, y0,tspan,args=(params,))
    Y_out[k]=yp[3,-1]+yp[4,-1]#max(yp[1,:])
#for j in range(0,p_len):
#    pyplot.figure(j)
#    pyplot.scatter(X_param[j,:], Y_out)
#this apparently ranks the data, could use rankdata from np.stat
Y_out_ranked=numpy.array(Y_out).argsort().argsort()
v_y_ranked= np.var(Y_out_ranked)
for j in range(0,p_len):
    X_param_ranked_all=numpy.array(X_param[j,:]).argsort().argsort()
for j in range(0,p_len):
    X_param_ranked=numpy.array(X_param[j,:]).argsort().argsort()
    v_pi_ranked=np.var(X_param_ranked, dtype=np.float64)
    cv_pi_y_ranked=np.cov(X_param_ranked, Y_out_ranked)
    SA_spearman[j]=cv_pi_y_ranked[0,1]/np.sqrt(v_pi_ranked*v_y_ranked)
#    SA_PRCC[j]=cv_pi_y_ranked
print(SA_spearman)
pyplot.figure(2) 
bar_width = 0.35 
pyplot.bar(np.arange(0,p_len)+bar_width/2.2,width=bar_width,height=SA_spearman, color='g')
pyplot.xticks(np.arange(0,p_len), param_list)

pyplot.figure(3) 
pyplot.scatter(X_param[0,:],Y_out)
pyplot.figure(4) 
pyplot.scatter(X_param[1,:],Y_out)
pyplot.figure(5) 
pyplot.scatter(X_param[2,:],Y_out)
pyplot.figure(6) 
pyplot.scatter(X_param[3,:],Y_out)
pyplot.figure(7) 
pyplot.scatter(X_param[7,:],Y_out)

