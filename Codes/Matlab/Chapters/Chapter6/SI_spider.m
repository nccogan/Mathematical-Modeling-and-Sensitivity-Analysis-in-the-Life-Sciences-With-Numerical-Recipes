%=========================================================
% Created on Thu Jun 16 07:30:47 2022
% Solves the SI model with vital dynamics. Uses Spider plot
% as the SI indice.
% @author: cogan
%=========================================================
clear;
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

params_names=( ['$r$', '$\kappa$', '$k$', '$\delta$', '$S_0$', '$I_0$'])
params=[.05, 100, .075, .3,99,1];
%===============================
% Plot one example
%===============================

[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [0 25],[params(5) params(6)]); 
plot(t, y_solution(:,1),'k','LineWidth',2)
hold on
plot(t, y_solution(:,2),'g','LineWidth',2)
legend(['Susceptible', 'Infected'], fontsize=16)
xlabel('Time', fontsize=16)
ylabel('Population', fontsize=16)
%savefig('SI_dyn.png')

%===============================
% Loop through the QoIs
%===============================

 params=[.05, 100, .075, .3,99,1]; %augment parameters for initial conditions

 Num_samples=5;
 QoI=zeros(2*Num_samples,6);
 for k  = 1:6
     for i = 1:2*Num_samples
        params(k)=.5*params(k)+(i-1)*.1*params(k);
        y0=[params(end-1),params(end)];
        [t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [0 500],[params(5) params(6)]); 
        QoI(i,k)=y_solution(end,1)/(y_solution(end,1)+y_solution(end,2));
        params=[.05, 100, .075, .3,99,1];
     end
 end
figure(2)
plot(QoI,'-o','Linewidth',2)
legend('$r$', '$\kappa$', '$k$', '$\delta$', '$S_0$', '$I_0$')
xticks(1:2*Num_samples+1)
xticklabels({'-50%', '-40%', '-30%', '-20%', '-10%', '0%', '10%', '20%', '30%', '40%', '50%'})
% pyplot.savefig('SI_spider.png')

%===============================
% Define the rhs functions
%===============================
function f1 = rhs_S(t,Y,params)
    S=Y(1);
    I=Y(2);
    r=params(1);
    kappa=params(2);
    k=params(3);
    delta=params(4);
    IC1=params(5);
    IC2 =params(6);
    f1 = -k*S*I+r*S*(kappa-S);
end


function f2 = rhs_I(t,Y,params)
    S=Y(1);
    I=Y(2);
    r=params(1);
    kappa=params(2);
    k=params(3);
    delta=params(4);
    IC1=params(5);
    IC2 =params(6);
    f2 = k*S*I-delta*I;
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
    g=[rhs_S(t,Y,params);rhs_I(t,Y,params)];
end

