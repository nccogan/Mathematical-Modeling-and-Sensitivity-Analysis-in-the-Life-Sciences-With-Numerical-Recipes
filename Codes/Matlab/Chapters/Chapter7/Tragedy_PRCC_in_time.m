% =============================================================================
% Created on Tue Oct 26 11:34:10 2021
% Cheater model. PRCC in time is implemented with pyDOE
% @author: cogan
% =============================================================================﻿
clear
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

Params_names=['$F$', '$S_0$', '$\gamma$', '$Y_1$', '$\mu_1$', '$K_1$', '$Y_2$', '$\mu_2$', '$K_2$',  '$\alpha$'];
% =============================================================================
% Note the differing initial conditions, here cheaters are present
% =============================================================================
params=[1, 1, 20,1, 5, .05, 1, 5, .05, .2];
tspan = linspace(0, 25, 500);
S_init=1;
P_init=0;
E_init=0.1;
B1_init=.2;
B2_init=0.02;
[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)], [S_init, P_init, E_init, B1_init, B2_init]); 
figure()

plot(t,y_solution(:,1), 'LineWidth',2)
hold on
plot(t,y_solution(:,2), 'LineWidth',2)
plot(t,y_solution(:,3), 'LineWidth',2)
plot(t,y_solution(:,4), 'LineWidth',2)
plot(t,y_solution(:,5), 'LineWidth',2)

xlabel('Time', fontsize=16)
ylabel('Concentration', fontsize=16)
legend('$S$', '$P$', '$E$', '$B_1$', '$B_2$', fontsize=16)
title('$\mu_1=5$ and $\bar{S}$ = 0.02 ' ,  'interpreter', 'latex',fontsize=16)

% =============================================================================
% Note the differing initial conditions, here cheaters are absent
% =============================================================================

params=[1, 1, 20,1, 5, .05, 1, 5, .05, .2];
tspan = linspace(0, 25, 500);
S_init=1;
P_init=0;
E_init=0.1;
B1_init=.2;
B2_init=0*0.02;
[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)], [S_init, P_init, E_init, B1_init, B2_init]); 
figure()

plot(t,y_solution(:,1), 'LineWidth',2)
hold on
plot(t,y_solution(:,2), 'LineWidth',2)
plot(t,y_solution(:,3), 'LineWidth',2)
plot(t,y_solution(:,4), 'LineWidth',2)
plot(t,y_solution(:,5), 'LineWidth',2)

xlabel('Time', fontsize=16)
ylabel('Concentration', fontsize=16)
legend('$S$', '$P$', '$E$', '$B_1$', '$B_2$', fontsize=16)
title('$\mu_1=5$ and $\bar{S}$ = 0 ' ,  'interpreter', 'latex',fontsize=16)

% =============================================================================
% Note the differing parameters, here mu_1 is small, no cheaters but the system 
% washout case is stable
% =============================================================================

params=[1, 1, 20,1, .5, .05, 1, 5, .05, .2];
tspan = linspace(0, 25, 500);
S_init=1;
P_init=0;
E_init=0.1;
B1_init=.2;
B2_init=0*0.02;
[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)], [S_init, P_init, E_init, B1_init, B2_init]); 
figure()

plot(t,y_solution(:,1), 'LineWidth',2)
hold on
plot(t,y_solution(:,2), 'LineWidth',2)
plot(t,y_solution(:,3), 'LineWidth',2)
plot(t,y_solution(:,4), 'LineWidth',2)
plot(t,y_solution(:,5), 'LineWidth',2)

xlabel('Time', fontsize=16)
ylabel('Concentration', fontsize=16)
legend('$S$', '$P$', '$E$', '$B_1$', '$B_2$', fontsize=16)
title('$\mu_1=.5$ and $\bar{S}$ = 0 ' ,  'interpreter', 'latex',fontsize=16)



% =============================================================================
% Note the differing parameters, here mu_1 is small, no cheaters but the system 
% washout case is almost stable
% =============================================================================

params=[1, 1, 20,1, 1.45, .05, 1, 5, .05, .2];
tspan = linspace(0, 25, 500);
S_init=1;
P_init=0;
E_init=0.1;
B1_init=.2;
B2_init=0*0.02;
[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)], [S_init, P_init, E_init, B1_init, B2_init]); 
figure()

plot(t,y_solution(:,1), 'LineWidth',2)
hold on
plot(t,y_solution(:,2), 'LineWidth',2)
plot(t,y_solution(:,3), 'LineWidth',2)
plot(t,y_solution(:,4), 'LineWidth',2)
plot(t,y_solution(:,5), 'LineWidth',2)

xlabel('Time', fontsize=16)
ylabel('Concentration', fontsize=16)
legend('$S$', '$P$', '$E$', '$B_1$', '$B_2$', fontsize=16)
title('$\mu_1=1.45$ and $\bar{S}$ = 0 ' ,  'interpreter', 'latex',fontsize=16)


% =============================================================================
%Partial correlation coefficient using partialcorr in Matlab
% =============================================================================

% =============================================================================
%This generates the LHS assuming normal distribution
%==============================================================================
Num_samples=100;
total_bacteria=zeros(1,Num_samples);
params=[1, 1, 20,1, 1.45, .05, 1, 5, .05, .2];
parameters=lhsnorm(params,diag(.05*params),Num_samples);

tspan = linspace(0, 25, 500);
S_init=1;
P_init=0;
E_init=0.1;
B1_init=.2;
B2_init=0.02;
PRCC_in_time=zeros([11,8]);
% =============================================================================
% Time loop for moving the ending time out time
% =============================================================================
tfinal_temp=[.001,.01,.1,1,5,10,30];
for k =1:8
    k
    tfinal=tfinal_temp(k);
    correlation_matrix=zeros([Num_samples,11]);
    for i = 1: Num_samples
        params_i=parameters(i,:);
        tspan = linspace(0, tfinal, 500);
        [t,y_solution]=ode45(@(t,Y) rhs(t,Y,params_i), [tspan(1),tspan(end)], [S_init, P_init, E_init, B1_init, B2_init]); 
        correlation_matrix(i,:)=[parameters(i,:),y_solution(end,4)];
    end

    
    PRCC=partialcorr(correlation_matrix);
    PRCC_in_time(:,k)=PRCC(:,end);
end
plot([.001,.01,.1,1,5,10,30],PRCC_in_time, 'LineWidth',2)
legend('$F$', '$S_0$', '$\gamma$', '$Y_1$', '$\mu_1$', '$K_1$', '$Y_2$', '$\mu_2$', '$K_2$',  '$\alpha$')
%==============================================
%Define the RHS of ODEs
%==============================================
function f1=rhs_S(t,Y,params)
    S=Y(1);
    P=Y(2);
    E=Y(3);
    B1=Y(4);
    B2=Y(5);
    F=params(1);
    S0=params(2);
    gamma=params(3);
    Yield1=params(4);
    mu1=params(5);
    K1=params(6);
    Yield2=params(7);
    mu2=params(8);
    K2=params(9);
    alpha=params(10);
     f1 = F*S0-F*S-gamma*S*E;
end

function f2=rhs_P(t,Y,params)
    S=Y(1);
    P=Y(2);
    E=Y(3);
    B1=Y(4);
    B2=Y(5);
    F=params(1);
    S0=params(2);
    gamma=params(3);
    Yield1=params(4);
    mu1=params(5);
    K1=params(6);
    Yield2=params(7);
    mu2=params(8);
    K2=params(9);
    alpha=params(10);
     f2 = gamma*S*E-1/Yield1*mu1*P/(K1+P)*B1-1/Yield2*mu2*P/(K2+P)*B2-F*P;
end

function f3=rhs_E(t,Y,params)
    S=Y(1);
    P=Y(2);
    E=Y(3);
    B1=Y(4);
    B2=Y(5);
    F=params(1);
    S0=params(2);
    gamma=params(3);
    Yield1=params(4);
    mu1=params(5);
    K1=params(6);
    Yield2=params(7);
    mu2=params(8);
    K2=params(9);
    alpha=params(10);
     f3 = alpha*mu1*P*B1/(K1+P)-F*E;
end

function f4=rhs_B1(t,Y,params)
    S=Y(1);
    P=Y(2);
    E=Y(3);
    B1=Y(4);
    B2=Y(5);
    F=params(1);
    S0=params(2);
    gamma=params(3);
    Yield1=params(4);
    mu1=params(5);
    K1=params(6);
    Yield2=params(7);
    mu2=params(8);
    K2=params(9);
    alpha=params(10);
     f4 = (1-alpha)*mu1*P*B1/(K1+P)-F*B1;
end

function f5=rhs_B2(t,Y,params)
    S=Y(1);
    P=Y(2);
    E=Y(3);
    B1=Y(4);
    B2=Y(5);
    F=params(1);
    S0=params(2);
    gamma=params(3);
    Yield1=params(4);
    mu1=params(5);
    K1=params(6);
    Yield2=params(7);
    mu2=params(8);
    K2=params(9);
    alpha=params(10);
     f5 =  mu2*P*B2/(K2+P)-F*B2;
end

function g=rhs(t,Y,params)
    S=Y(1);
    P=Y(2);
    E=Y(3);
    B1=Y(4);
    B2=Y(5);
    F=params(1);
    S0=params(2);
    gamma=params(3);
    Yield1=params(4);
    mu1=params(5);
    K1=params(6);
    Yield2=params(7);
    mu2=params(8);
    K2=params(9);
    alpha=params(10);
    g = [rhs_S(t,Y,params);rhs_P(t,Y,params);rhs_E(t,Y,params);rhs_B1(t,Y,params);rhs_B2(t,Y,params)];
end


