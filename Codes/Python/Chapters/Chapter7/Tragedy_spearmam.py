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
We use the PRCC version in the package pingouin. We also show native commands for the 
ranked correlation coefficienct to compare the estimates when discounting. 

This can be compared to the code in cheater.py that uses the package pyDOE
@author: cogan
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import numpy.matlib
from scipy.integrate import odeint
import matplotlib.pyplot as pyplot
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.integrate import solve_ivp
import scipy as scipy
pyplot.close('all')
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible

#Define the RHS of ODEs
def f1(t,Y,params):
    S,P,E,X1,X2=Y
    S0,D,k1,k2,Yield,mu,K,q=params
    return D*(S0-S)-(k1*S*E-k2*P)
def f2(t,Y,params):
    S,P,E,X1,X2=Y
    S0,D,k1,k2,Yield,mu,K,q=params
    return (k1*S*E-k2*P)-1/Yield*(X1+X2)*mu*P/(K+P)-D*P
def f3(t,Y,params):
    S,P,E,X1,X2=Y
    S0,D,k1,k2,Yield,mu,K,q=params
    return (1-q)*X1*mu*P/(K+P)-D*E
def f4(t,Y,params):
    S,P,E,X1,X2=Y
    S0,D,k1,k2,Yield,mu,K,q=params
    return X1*q*mu*P/(K+P)-D*X1
def f5(t,Y,params):
    S,P,E,X1,X2=Y
    S0,D,k1,k2,Yield,mu,K,q=params
    return X2*mu*P/(K+P)-D*X2
def rhs(t,Y,params):
    S,P,E,X1,X2=Y
    S0,D,k1,k2,Yield,mu,K,q=params
    return [f1(t,Y,params),f2(t,Y,params),f3(t,Y,params),f4(t,Y,params),f5(t,Y,params)]
#####################from scipy.integrate import ode
param_list=['S','D','k1','k2','Yield','mu','K','q']
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
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Concentration', fontsize=16)
pyplot.legend(['$S$', '$P$', '$E$', '$X_1$', '$X_2$' ], fontsize=16)

####Make the parameter matrix
N=500 # number of samples
del_p=.1 # parameter variation e.g. 10%
p_len=len(params_0)
X_param=np.zeros((p_len,N))
Y_out=np.zeros((N))
SA_pearson=np.zeros(p_len)
SA_spearman=np.zeros(p_len)
SA_PRCC=np.zeros(p_len)

tspan = np.linspace(0, 2, 500)

for i in range(0,p_len):
    #one parameter
    X_param[i,:]=np.random.uniform(params_0[i]*(1-del_p),params_0[i]*(1+del_p),N) 

for k in range(0,N):
    params=X_param[:,k]  
    yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan1[0],tspan1[-1]], y0, method='RK45',t_eval=tspan1)
    Y_out[k]=max(yp.y[1,:])
#    Y_out[k]=max(yp[1,:])
v_y=np.var(Y_out, dtype=np.float64)


Y_out_ranked=numpy.array(Y_out).argsort().argsort()
v_y_ranked= np.var(Y_out_ranked)
for j in range(0,p_len):
    X_param_ranked_all=numpy.array(X_param[j,:]).argsort().argsort()
for j in range(0,p_len):
    X_param_ranked=numpy.array(X_param[j,:]).argsort().argsort()
    v_pi_ranked=np.var(X_param_ranked, dtype=np.float64)
    cv_pi_y_ranked=np.cov(X_param_ranked, Y_out_ranked)
    SA_spearman[j]=cv_pi_y_ranked[0,1]/np.sqrt(v_pi_ranked*v_y_ranked)
print(SA_spearman)
pyplot.figure(2) 
bar_width = 0.35 
pyplot.bar(np.arange(0,p_len)+bar_width/2.2,width=bar_width,height=SA_spearman, color='g')
pyplot.xticks(np.arange(0,p_len), param_list)




# =============================================================================
#Alternative implementation
# #pip install pingouin
# =============================================================================
import pandas as pd
v=['S','D','k1','k2','Yield','mu','K','q','X'];
df = pd.DataFrame(np.vstack((X_param, Y_out.T)).T, columns = ['S','D','k1','k2','Yield','mu','K','q','X'])
prcc=np.zeros(p_len)
import pingouin as pg
for i in range(0,p_len):
    list=v[0:i]+v[i+1:-2]
    prcc[i]=pg.partial_corr(data=df, x=v[i], y='X', covar=list,method='spearman').round(3).r

pyplot.bar(np.arange(0,p_len)-bar_width/2.2,width=bar_width,height=SA_spearman, color='g')

pyplot.bar(np.arange(0,p_len)+bar_width/2.2,width=bar_width,height=prcc, color='b')


# =============================================================================
# Scatter plots to check for monotonicity
# =============================================================================
Y_out_0=np.empty(N)
for k in range(0,N):
    params=[ X_param[0,k],1, 20, .005, 1, 5,.05,.8]
    yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan1[0],tspan1[-1]], y0, method='RK45',t_eval=tspan1)
    Y_out_0[k]=max(yp.y[1,:])
pyplot.figure(3) 

pyplot.scatter(X_param[0,:],Y_out_0)
pyplot.xlabel('S', fontsize=16)
pyplot.ylabel('QoI', fontsize=16)

Y_out_1=np.empty(N)
for k in range(0,N):
    params=[ 1, X_param[1,k], 20, .005, 1, 5,.05,.8]
    yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan1[0],tspan1[-1]], y0, method='RK45',t_eval=tspan1)
    Y_out_1[k]=max(yp.y[1,:])


pyplot.figure(4) 
pyplot.scatter(X_param[1,:],Y_out_1)
pyplot.xlabel('D', fontsize=16)
pyplot.ylabel('QoI', fontsize=16)

Y_out_2=np.empty(N)
for k in range(0,N):
    params=[ 1, 1, X_param[2,k], .005, 1, 5,.05,.8]
    yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan1[0],tspan1[-1]], y0, method='RK45',t_eval=tspan1)
    Y_out_2[k]=max(yp.y[1,:])


pyplot.figure(5) 
pyplot.scatter(X_param[2,:],Y_out_2)
pyplot.xlabel('k_1', fontsize=16)
pyplot.ylabel('QoI', fontsize=16)


Y_out_3=np.empty(N)
for k in range(0,N):
    params=[ 1, 1, 20, X_param[3,k], 1, 5,.05,.8]
    yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan1[0],tspan1[-1]], y0, method='RK45',t_eval=tspan1)
    Y_out_3[k]=max(yp.y[1,:])

pyplot.figure(6) 
pyplot.scatter(X_param[3,:],Y_out_3)
pyplot.xlabel('k_2', fontsize=16)
pyplot.ylabel('QoI', fontsize=16)



Y_out_4=np.empty(N)
for k in range(0,N):
    params=[ 1, 1, 20, .005, X_param[4,k], 5,.05,.8]
    yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan1[0],tspan1[-1]], y0, method='RK45',t_eval=tspan1)
    Y_out_4[k]=max(yp.y[1,:])


pyplot.figure(7) 
pyplot.scatter(X_param[4,:],Y_out_4)
pyplot.xlabel('Yield', fontsize=16)
pyplot.ylabel('QoI', fontsize=16)


Y_out_5=np.empty(N)
for k in range(0,N):
    params=[ 1, 1, 20, .005, 1, X_param[5,k],.05,.8]
    yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan1[0],tspan1[-1]], y0, method='RK45',t_eval=tspan1)
    Y_out_5[k]=max(yp.y[1,:])

pyplot.figure(8) 
pyplot.scatter(X_param[5,:],Y_out_5)
pyplot.xlabel('$/mu$', fontsize=16)
pyplot.ylabel('QoI', fontsize=16)


Y_out_6=np.empty(N)
for k in range(0,N):
    params=[ 1, 1, 20, .005, 1, 5,X_param[6,k],.8]
    yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan1[0],tspan1[-1]], y0, method='RK45',t_eval=tspan1)
    Y_out_6[k]=max(yp.y[1,:])

pyplot.figure(9) 
pyplot.scatter(X_param[6,:],Y_out_6)
pyplot.xlabel('K', fontsize=16)
pyplot.ylabel('QoI', fontsize=16)

Y_out_7=np.empty(N)
for k in range(0,N):
    params=[ 1, 1, 20, .005, 1, 5,.05,X_param[7,k]]
    yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan1[0],tspan1[-1]], y0, method='RK45',t_eval=tspan1)
    Y_out_6[k]=max(yp.y[1,:])

pyplot.figure(10) 
pyplot.scatter(X_param[7,:],Y_out_6)
pyplot.xlabel('q', fontsize=16)
pyplot.ylabel('QoI', fontsize=16)


