
% Created on Mon Jan  3 17:42:39 2022
% Model of SIR
% @author: Cogan

clear;
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');



p_names=[ 'k','gamma'];
params=[ 0.005, 1];
N=1000;
tspan = linspace(0, 10, 500);

[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)], [N,1,0]); 
plot(t, y_solution(:,1),'k','LineWidth',2)
hold on
plot(t, y_solution(:,2),'k','LineWidth',2)
plot(t, y_solution(:,3),'k','LineWidth',2)

legend(['Susceptible','Infected','Recovered'], fontsize=16)
xlabel('Time', fontsize=16)
ylabel('Population', fontsize=16)

figure(2)
params1=[ 0.005,2];

[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)], [N,1,0]); 
plot(t, y_solution(:,1),'k','LineWidth',2)
hold on
plot(t, y_solution(:,2),'b','LineWidth',2)
plot(t, y_solution(:,3),'c','LineWidth',2)

legend(['Susceptible','Infected','Recovered'], fontsize=16)
xlabel('Time', fontsize=16)
ylabel('Population', fontsize=16)


%Generate alpha plot
figure(3)
params1=[ 0.01, 0.5];
[t2,y_solution2]=ode45(@(t,Y) rhs(t,Y,params1), [tspan(1),tspan(end)], [.75*N,1,0]); 
[t3,y_solution3]=ode45(@(t,Y) rhs(t,Y,params1), [tspan(1),tspan(end)], [.5*N,1,0]); 
[t4,y_solution4]=ode45(@(t,Y) rhs(t,Y,params1), [tspan(1),tspan(end)], [.25*N,1,0]); 
[t5,y_solution5]=ode45(@(t,Y) rhs(t,Y,params1), [tspan(1),tspan(end)], [.05*N,1,0]); 

plot(t2, y_solution2(:,2),'k','LineWidth',2)
hold on
plot(t3, y_solution3(:,2),'b','LineWidth',2)
plot(t4, y_solution4(:,2),'c','LineWidth',2)
plot(t5, y_solution5(:,2),'m','LineWidth',2)
legend('$\alpha = .75$','$\alpha = .5$','$\alpha = .25$','$\alpha = .05$', fontsize=16)
xlabel('Time', fontsize=16)
ylabel('Population', fontsize=16)


%===============================
%      Tornado plot
% we need low and high values of 
% parameters with results Qoi_low 
% and QoI_high
%===============================




params_t=[.1, 2000, .01, .01, .5, 1000-10,10,0]; %augment parameters for initial conditions
params_min_t=[.1, 2000, .01, .01, .5, 1000-10,10,0];
params_max_t=[.1, 2000, .01, .01, .5, 1000-10,10,0];
QoI_min=zeros(1,7);
QoI_max=zeros(1,7);

for k =1:7
    params_min_t(k)=.5*params_min_t(k);
    params_max_t(k)=1.5*params_max_t(k);
    [t2_min_t,y_solution_min_t]=ode45(@(t,Y) rhs_t(t,Y,params_min_t), [tspan(1),tspan(end)], [params_t(end-2),params_t(end-1), params_t(end)]); 
    [t2_max_t,y_solution_max_t]=ode45(@(t,Y) rhs_t(t,Y,params_max_t), [tspan(1),tspan(end)], [params_t(end-2),params_t(end-1), params_t(end)]); 
    QoI_min(k)=y_solution_min_t(end,1)/(y_solution_min_t(end,1)+y_solution_min_t(end,2));
    QoI_max(k)=y_solution_max_t(end,1)/(y_solution_max_t(end,1)+y_solution_max_t(end,2));
    params_t=[.1, 2000, .01, .01, .5, 1000-10,10,0]; %augment parameters for initial conditions
    params_min_t=[.1, 2000, .01, .01, .5, 1000-10,10,0];
    params_max_t=[.1, 2000, .01, .01, .5, 1000-10,10,0];
end

%========================
%Make the tornado plots
%========================
[QoI_min,ind]=sort(QoI_min,'descend');
QoI_max=QoI_max(ind);

figure
h = barh(QoI_min);
hold on
barh(-QoI_max,'r')

xticklabels({[[max(QoI_max),max(QoI_max)/2, 0,max(QoI_max)/2, max(QoI_max)]]})
yticklabels({'$r$','$\kappa$','$\gamma$','$\delta$','$S_0$','$I_0$','$R_0$'})

%=========================
%Define the RHS of ODEs
%=========================
function f1 = rhs_S(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    k=params(1);
    gamma=params(2);
    f1 = -k*S*I;
end

function f2 = rhs_I(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    k=params(1);
    gamma=params(2);
    f2 = k*S*I-gamma*I;
end

function f3 = rhs_R(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    k=params(1);
    gamma=params(2);
    f3 = gamma*I;
end
function g = rhs(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    k=params(1);
    gamma=params(2);
    g = [rhs_S(t,Y,params);rhs_I(t,Y,params);rhs_R(t,Y,params)];
end

%=========================
%Define the RHS of ODEs
%for tornado plots
%=========================
function f1 = rhs_S_t(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    r=params(1);
    kappa=params(2);
    k=params(3);
    gamma=params(4);
    delta=params(5);
    IC_S=params(6);
    IC_I=params(7);
    IC_R=params(8);
    f1 = -k*S*I;
end

function f2 = rhs_I_t(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    r=params(1);
    kappa=params(2);
    k=params(3);
    gamma=params(4);
    delta=params(5);
    IC_S=params(6);
    IC_I=params(7);
    IC_R=params(8);
    f2 = k*S*I-gamma*I;
end

function f3 = rhs_R_t(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    r=params(1);
    kappa=params(2);
    k=params(3);
    gamma=params(4);
    delta=params(5);
    IC_S=params(6);
    IC_I=params(7);
    IC_R=params(8);
    f3 = gamma*I;
end
function g = rhs_t(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    r=params(1);
    kappa=params(2);
    k=params(3);
    gamma=params(4);
    delta=params(5);
    IC_S=params(6);
    IC_I=params(7);
    IC_R=params(8);
    g = [rhs_S_t(t,Y,params);rhs_I_t(t,Y,params);rhs_R_t(t,Y,params)];
end













