% =============================================================================﻿
% Created on Tue Oct 26 11:34:10 2021
% Example of freezing parameters, Freter model. The idea is to run the model
% will all parameters varying and comparing the statistics of the output to 
% samples of the model with only the sensitive parameters varying.
% 
% This is one standard use of sensitivity -- it can reduce the parameter space needed
% to provide accurate predictions.
% @author: cogan
% =============================================================================

clear
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');





% =============================================================================
%This generates the LHS assuming normal distribution
%We use PRCC to determine the sensitivities
%==============================================================================
Params_names=['$N_0$', '$Y$', '$\mu$', '$K_n$', '$F$', '$B_0$', '$\alpha$', '$k_\alpha$', '$V$', '$A$', '$\beta$', '$K_b$'] ;
params=[1, .1, 1, .5, 1, .1, .1, .1, 10, 1, .2, .5];
tspan = linspace(0, 40, 500);
N_init=0;
Bu_init=.1;
Bb_init=0;
Num_samples=500;
total_bacteria=zeros(1,Num_samples);
parameters=abs(lhsnorm(params,diag(.05*params),Num_samples));

for i =1:Num_samples
    params_i=parameters(i,:);
    [t,y_solution]=ode45(@(t,Y) rhs(t,Y,params_i), [tspan(1),tspan(end)], [N_init, Bu_init, Bb_init]); 
    total_bacteria(i)=trapz(t,y_solution(:,2))/trapz(t,y_solution(:,3)); %trapz is a relatively standard  implementation of the trapezoidal rule for integration
end
cc=zeros(1,12);
cc1=zeros(1,12);

for j = 1:length(params)
   cc_mat=corrcoef(parameters(:,j), total_bacteria);
   cc(j)=cc_mat(1,end);
end
    
for j = 1:length(params)
    cc1_mat(j)=corr(parameters(:,j), total_bacteria','Type','Spearman');
    cc1(j)=cc1_mat(1,end);
end


figure()
bar_width=.35;
p_names=['$N_0$', '$F$', '$Y$', '$\mu$', '$K_n$' ];
hold on
bar([1:length(params)]-.5*bar_width,cc1,bar_width ,'b')
bar([1:length(params)]+.5*bar_width,cc, bar_width,'k')
xticks(1:length(params))
xticklabels({'$N_0$', '$Y$', '$\mu$', '$K_n$', '$F$', '$B_0$', '$\alpha$', '$k_\alpha$', '$V$', '$A$', '$\beta$', '$K_b$'})
legend('Pearson', 'Spearman')

% =============================================================================
% It is relatively simple to see that several of these are not sensitive
% we freeze those by keeping them at their nominal value (line 111 below)
% =============================================================================
% =============================================================================
% Frozen histograms
% =============================================================================
total_bacteria_frozen=zeros(1,Num_samples);

for i =1:Num_samples
    params_i=parameters(i,:);
    params_i(1:6)=params(1:6);
    [t,y_solution]=ode45(@(t,Y) rhs(t,Y,params_i), [tspan(1),tspan(end)], [N_init, Bu_init, Bb_init]); 
    total_bacteria_frozen(i)=trapz(t,y_solution(:,2))/trapz(t,y_solution(:,3)); %trapz is a relatively standard  implementation of the trapezoidal rule for integration
end

figure()

histogram(total_bacteria(find(total_bacteria<=10&total_bacteria>.01)),50) %discounting some large values
hold on
histogram(total_bacteria_frozen(find(total_bacteria_frozen<=10&total_bacteria_frozen>.01)),50) %discounting some large values

legend('Full', 'Frozen',fontsize=16)
xlabel('QoI',fontsize=16)
ylabel('Frequency',fontsize=16)



% =============================================================================
% Define the RHS of ODEs
% =============================================================================
function f1=rhs_N(t,Y,params)
    N=Y(1);
    Bu=Y(2);
    Bb=Y(3);
    N0=params(1);
    Yield=params(2);
    mu=params(3);
    K_n=params(4);
    F=params(5);
    B0=params(6);
    alpha=params(7);
    K_alpha=params(8);
    V=params(9);
    A=params(10);
    beta=params(11);
    K_b=params(12);
    f1= N0*F-1/Yield*mu*N/(K_n+N)*(Bu+Bb)-F*N;
end

function f2=rhs_Bu(t,Y,params)
    N=Y(1);
    Bu=Y(2);
    Bb=Y(3);
    N0=params(1);
    Yield=params(2);
    mu=params(3);
    K_n=params(4);
    F=params(5);
    B0=params(6);
    alpha=params(7);
    K_alpha=params(8);
    V=params(9);
    A=params(10);
    beta=params(11);
    K_b=params(12);
    f2 =  B0*F-alpha/(K_alpha+Bb)*Bu+V/A*beta*Bb ...
            -F*Bu +(1-Bb/(K_b+Bb))*1/Yield*mu*N/(K_n+N)*Bb-F*Bu;
end

function f3=rhs_Bb(t,Y,params)
    N=Y(1);
    Bu=Y(2);
    Bb=Y(3);
    N0=params(1);
    Yield=params(2);
    mu=params(3);
    K_n=params(4);
    F=params(5);
    B0=params(6);
    alpha=params(7);
    K_alpha=params(8);
    V=params(9);
    A=params(10);
    beta=params(11);
    K_b=params(12);
    f3 =  A/V*alpha/(K_alpha+Bb)*Bu-beta*Bb ...
        +A/V*Bb/(K_b+Bb)*1/Yield*mu*N/(K_n+N)*Bb;
end

function g=rhs(t,Y,params)
    N=Y(1);
    Bu=Y(2);
    Bb=Y(3);
    N0=params(1);
    Yield=params(2);
    mu=params(3);
    K_n=params(4);
    F=params(5);
    B0=params(6);
    alpha=params(7);
    K_alpha=params(8);
    V=params(9);
    A=params(10);
    beta=params(11);
    K_b=params(12);
    g = [rhs_N(t,Y,params);rhs_Bu(t,Y,params);rhs_Bb(t,Y,params)];
end

