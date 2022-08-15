"""
Created on Mon Jan  3 17:42:39 2022
Model of SIR
@author: Cogan
"""


from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.integrate import odeint
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
import random
pyplot.close('all')
pyplot.style.use('tableau-colorblind10')  #Fixes colorscheme to be accessible



#Define the RHS of ODEs
def f1(t,Y,params):
    S,I,R=Y
    k, gamma = params
    return -k*S*I
def f2(t,Y,params):
    S,I,R=Y
    k, gamma = params
    return k*S*I-gamma*I
def f3(t,Y,params):
    S,I,R=Y
    k, gamma = params
    return gamma*I
def rhs(t,Y,params):
    S,I,R=Y
    k, gamma = params
    return [f1(t,Y,params),f2(t,Y,params),f3(t,Y,params)]
# same as usual: define f and g, rhs returns the right hand side
# parameters are from the paper, beware the order though
p_names=( 'k','gamma')
params=np.array([ 0.005, 1])
N=1000
tspan = np.linspace(0, 10, 500)
yp= solve_ivp(lambda t,Y: rhs(t,Y,params), [tspan[0],tspan[-1]], [N,1,0], method='RK45',t_eval=tspan)
pyplot.figure(1)

pyplot.plot(tspan,yp.y[0])
pyplot.plot(tspan,yp.y[1])
pyplot.plot(tspan,yp.y[2])
pyplot.legend(['Susceptible','Infected','Recovered'], fontsize=16)
pyplot.legend(['Susceptible','Infected','Recovered'], fontsize=16)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Population', fontsize=16)
#pyplot.savefig('SIR_high.png')

pyplot.figure(2)
params1=np.array([ 0.005,2])
yp1= solve_ivp(lambda t,Y: rhs(t,Y,params1), [tspan[0],tspan[-1]], [N,1,0], method='LSODA',t_eval=tspan)
pyplot.plot(tspan,yp1.y[0])
pyplot.plot(tspan,yp1.y[1])
pyplot.plot(tspan,yp.y[2])
pyplot.legend(['Susceptible','Infected','Recovered'], fontsize=16)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Population', fontsize=16)
#pyplot.savefig('SIR_low.png')


###GENERATE alpha plot
pyplot.figure(3)
params1=np.array([ 0.01, 0.5])
yp2= solve_ivp(lambda t,Y: rhs(t,Y,params1), [tspan[0],tspan[-1]], [.75*N,1,0], method='RK45',t_eval=tspan)
yp3= solve_ivp(lambda t,Y: rhs(t,Y,params1), [tspan[0],tspan[-1]], [.5*N,1,0], method='RK45',t_eval=tspan)
yp4= solve_ivp(lambda t,Y: rhs(t,Y,params1), [tspan[0],tspan[-1]], [.25*N,1,0], method='RK45',t_eval=tspan)
yp5= solve_ivp(lambda t,Y: rhs(t,Y,params1), [tspan[0],tspan[-1]], [.05*N,1,0], method='RK45',t_eval=tspan)

pyplot.plot(tspan,yp2.y[1])
pyplot.plot(tspan,yp3.y[1])
pyplot.plot(tspan,yp4.y[1])
pyplot.plot(tspan,yp5.y[1])

pyplot.legend([r'$\alpha = .75$',r'$\alpha = .5$',r'$\alpha = .25$',r'$\alpha = .05$'], fontsize=16)
pyplot.xlabel('Time', fontsize=16)
pyplot.ylabel('Population', fontsize=16)
#pyplot.savefig('SIR_alpha.png')


###############################################################################
#Tornado plot
#we need low and high values of parameters with results Qoi_low and QoI_high
###############################################################################




#Define the RHS of ODEs for the tornado plots. We have expanded the number of parameters
def f1_t(t,Y,params):
    S,I,R=Y
    r, kappa, k, gamma, delta, IC_S, IC_I, IC_R = params
    return r*S*(kappa-S)-k*S*I
def f2_t(t,Y,params):
    S,I,R=Y
    r, kappa, k, gamma, delta , IC_S, IC_I, IC_R= params
    return k*S*I-gamma*I-delta*I
def f3_t(t,Y,params):
    S,I,R=Y
    r, kappa, k, gamma, delta , IC_S, IC_I, IC_R= params
    return gamma*I
def rhs_t(t,Y,params):
    S,I,R=Y
    r, kappa, k, gamma, delta , IC_S, IC_I, IC_R= params
    return [f1_t(t,Y,params),f2_t(t,Y,params),f3_t(t,Y,params)]


tspan = np.linspace(0, 50, 1000)

params_t=[.1, 2000, .01, .01, .5, 1000-10,10,0] #augment parameters for initial conditions
params_min_t=[.1, 2000, .01, .01, .5, 1000-10,10,0]
params_max_t=[.1, 2000, .01, .01, .5, 1000-10,10,0]
QoI_min=np.zeros([7])
QoI_max=np.zeros([7])

for k in np.arange(0,7):
    params_min_t[k]=.5*params_min_t[k]
    params_max_t[k]=1.5*params_max_t[k]
    y_solution_min_t = solve_ivp(lambda t,Y: rhs_t(t,Y,params_min_t),  [tspan[0],tspan[-1]],[params_t[-3],params_t[-2], params_t[-1]],method='RK45',t_eval=tspan)
    y_solution_max_t = solve_ivp(lambda t,Y: rhs_t(t,Y,params_max_t),  [tspan[0],tspan[-1]],[params_t[-3],params_t[-2], params_t[-1]],method='RK45',t_eval=tspan)     
    QoI_min[k]=y_solution_min_t.y[0][-1]/(y_solution_min_t.y[0][-1]+y_solution_min_t.y[1][-1])
    QoI_max[k]=y_solution_max_t.y[0][-1]/(y_solution_max_t.y[0][-1]+y_solution_max_t.y[1][-1])
    params_t=[.1, 2000, .01, .01, .5, 1000-10,10,0] #augment parameters for initial conditions
    params_min_t=[.1, 2000, .01, .01, .5, 1000-10,10,0]
    params_max_t=[.1, 2000, .01, .01, .5, 1000-10,10,0]


#Make the tornado plots
pos_label = (['$r$','$\kappa$','$\gamma$','$\delta$','$S_0$','$I_0$','$R_0$'])#
pos=np.arange(7) + .5

fig, (ax_left, ax_right) = pyplot.subplots(ncols=2)
ax_left.barh(pos, QoI_min, align='center')
ax_left.set_yticks([])
ax_left.invert_xaxis()
ax_right.barh(pos, QoI_max, align='center')
ax_right.set_yticks(pos)
ax_right.set_yticklabels(pos_label, ha='center', x=-1)
ax_right.set_xlim((0,.0075))
ax_left.set_xlim((0.0075,0))
pyplot.savefig('SIR_tornado.png')























