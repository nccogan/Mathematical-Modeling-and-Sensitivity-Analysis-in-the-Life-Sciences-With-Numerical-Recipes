%====================================﻿
%Created on Tue Oct 26 11:34:10 2021
%SIRS model with cobwebbing
%@author: cogan
%====================================

clear
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

p_names=['$k$','$\alpha$', '$\gamma$', '$S_0$','$I_0$','$R_0$' ];

%With no reversion from recovered to susceptible  (alpha=0) there is only one wave
params=[ 0.001, 0, .2, 500, 10, 0];
tspan = linspace(0, 100, 500);
[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)], [params(end-2),params(end-1),params(end)]); 
figure(1)

plot(t,y_solution(:,1), 'LineWidth',2)
hold on
plot(t,y_solution(:,2), 'LineWidth',2)
plot(t,y_solution(:,3), 'LineWidth',2)

legend(['Susceptible','Infected','Recovered'], fontsize=16)
xlabel('Time', fontsize=16)
ylabel('Population', fontsize=16)


%Nonzero alpha leades to oscillating solutions
params=[ 0.001, 0.0125, .2, 500, 10, 0];
tspan = linspace(0, 200, 500);
[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)], [params(end-2),params(end-1),params(end)]); 
figure(2)

plot(t,y_solution(:,1), 'LineWidth',2)
hold on
plot(t,y_solution(:,2), 'LineWidth',2)
plot(t,y_solution(:,3), 'LineWidth',2)

legend(['Susceptible','Infected','Recovered'], fontsize=16)
xlabel('Time', fontsize=16)
ylabel('Population', fontsize=16)

%With a faster rate of waning antigens, the disease rapidly becomes endemic
params=[ 0.001, 0.25, .2, 500, 10, 0];
tspan = linspace(0, 200, 500);
[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)], [params(end-2),params(end-1),params(end)]); 
figure(3)

plot(t,y_solution(:,1), 'LineWidth',2)
hold on
plot(t,y_solution(:,2), 'LineWidth',2)
plot(t,y_solution(:,3), 'LineWidth',2)

legend(['Susceptible','Infected','Recovered'], fontsize=16)
xlabel('Time', fontsize=16)
ylabel('Population', fontsize=16)


% =============================================================================
% Set-up cobwebbing 
% =============================================================================

% =============================================================================
% set-up samples
% =============================================================================
params=[ 0.001, 0.01, .2, 500, 10, 1];
tspan = linspace(0, 50, 500);
Num_samples=200;
param_level=linspace(-1, 1, Num_samples); %this iterates the scaling of the parameter samples
parameter_list=zeros([Num_samples,6]);
% =============================================================================
% Looop through each parameter and each value -- arrange in a matrix, use this
% to generate samples.  Note for loops are pretty fast in 
% python, but typically not in Matlab. We could make this more efficient. 
% =============================================================================
for k = 1:6
    for i = 1: Num_samples
        parameter_list(i,k)=(1+randsample(param_level,1))*params(k);
    end
end
        
for i = 1: Num_samples    
    [t,y_solution]=ode45(@(t,Y) rhs(t,Y,parameter_list(i,1:6)), [tspan(1),tspan(end)], [params(end-2),params(end-1),params(end)]); 
    parameter_list(i,end)=y_solution(end,2);
end

% =============================================================================
% Find all QoI (here it is the infected population at time 200) larger than a threshold
% (here 150) and sort samples into large and small infected populations
% =============================================================================
Large=find(parameter_list(:,end)>150); %find all QoI>150
Small=find(parameter_list(:,end)<1); %find all QoI<1

% =============================================================================
% Scale parameters between -1 and 1
% =============================================================================
parameter_list_scaled=zeros(Num_samples,6);
scaled=zeros(1,6);

for i  = 1: Num_samples
    scaled=parameter_list(i,1:6)./[params(1:5),1];
    parameter_list_scaled(i,1:6)=scaled-1;
end

    
% =============================================================================
%   Generate Spider plot
% =============================================================================
figure(4)
hold on

for i = 1 : length(Large)
     plot(parameter_list_scaled(Large(i),1:end-1),'co-', 'LineWidth',2)
end


for i = 1 : length(Small)
    plot(parameter_list_scaled(Small(i),1:end-1),'k*-', 'LineWidth',2)
end
h=plot(parameter_list_scaled(Large(end),1:end-1),'co-', 'LineWidth',2);

g=plot(parameter_list_scaled(Small(end),1:end-1),'k*-', 'LineWidth',2);

xticks(1:5)
xticklabels({'$k$','$\alpha$', '$\gamma$', '$S_0$','$I_0$','$R_0$'})
legend([h,g],'High', 'Low')
% pyplot.savefig('SIRS_cobweb.png')         
% 
% 
% 

% =============================================================================
%     Define right-hand-sides
% =============================================================================
function f1 = rhs_S(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    k=params(1);
    alpha=params(2);
    gamma=params(3);
    S0=params(4);
    I0=params(5);
    R0=params(6);
    f1 = -k*S*I+alpha*R;
end

function f2 = rhs_I(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    k=params(1);
    alpha=params(2);
    gamma=params(3);
    S0=params(4);
    I0=params(5);
    R0=params(6);
    f2 = k*S*I-gamma*I;
end

function f3 = rhs_R(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    k=params(1);
    alpha=params(2);
    gamma=params(3);
    S0=params(4);
    I0=params(5);
    R0=params(6);
    f3 = gamma*I-alpha*R;
end

function g = rhs(t,Y,params)
    S=Y(1);
    I=Y(2);
    r=params(1);
    kappa=params(2);
    k=params(3);
    delta=params(4);
    IC1=params(5);
    IC2 =params(6);
    g=[rhs_S(t,Y,params);rhs_I(t,Y,params);rhs_R(t,Y,params)];
end

