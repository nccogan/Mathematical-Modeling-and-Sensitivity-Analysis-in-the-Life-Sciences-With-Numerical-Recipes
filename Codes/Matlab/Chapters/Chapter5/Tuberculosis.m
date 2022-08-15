%============================================﻿
% Created on Fri Dec 31 22:28:06 2021
% TB Model
% @author: cogan
% B= Free bacteria
% Q = Dormant Bacteria
% X= Immune response
% r,g,h,f,f,g,a,s,k,d=params
%=============================================
clear;
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

tstart=0;
tstop=160;
tspan = linspace(tstart, tstop, 200);
y0 = [1,1,1];
params=[1,.001,.5,.1,.1,1,100,.1]; 

[t,y_simulation]=ode45(@(t,Y) rhs(t,Y,params), [tstart tstop],y0); 
figure(1)
plot(t,log(y_simulation(:,1)),'b',LineWidth=2)
hold on
plot(t,log(y_simulation(:,1)+y_simulation(:,2)),'m',LineWidth=2)
plot(t,log(y_simulation(:,3)),'k',LineWidth=2)
legend('B','B+Q','X',fontsize=16)
xlabel('Time',fontsize=16)
ylabel('log(density)',fontsize=16)

% %=====================
% %function to calculate 
% % %the time of maximum
% %=====================
Num_samples=1000;
params_max=1.1*params;
params_min=.9*params;
time_max=zeros(1,Num_samples);
 for s =1:Num_samples  
     for k = 1:length(params) 
         params(k) = params_min(k) + (params_max(k)-params_min(k)).*rand(1,1);
        [t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tstart tstop],y0); 
        [val, idx] = max(y_solution(:,1));
        time_max(s)=t(idx);
     end
 end

figure(2)
hist(time_max,15)
xlabel('Time of Maximum',fontsize=16)
ylabel('Frequency',fontsize=16)
% pyplot.savefig('time_max_hist.png')
 
%=====================
%Local Sensitivity
%=====================
Num_samples=10;
params_max=1.2*params;
params_min=.9*params;
time_max=zeros(1,Num_samples);
s=zeros(1,8); %place holder for der. of QoI
r=zeros(1,8); %place holder for params

time_max_min=zeros(1,8);
time_max_max=zeros(1,8);
params_small=params;
params_large=params;
 
for k = 1:8 
    params_small(k)=params_min(k);
    [t,y_solution_small]=ode45(@(t,Y) rhs(t,Y,params_small), [tstart tstop],y0); 
    [val_min, idx_min] = max(y_solution_small(:,1));
    time_max_min(k)=t(idx_min);

    params_large(k)=params_max(k);
    [t,y_solution_large]=ode45(@(t,Y) rhs(t,Y,params_large), [tstart tstop],y0); 
    [val_max, idx_max] = max(y_solution_large(:,1));
    time_max_max(k)=t(idx_max);

    s(k)=(time_max_max(k)-time_max_min(k))/time_max_max(k);
    r(k)=(params_large(k)-params_small(k))/params(k);
    params_small=params;
    params_large=params;   
end
figure(3)

bar(1:8,s./r)
xticks(1:8)
xticklabels({'$r$', '$g$', '$h$', '$f$', '$a$','$s$','$k$','$d$'})


%=====================
%Define the functions
%=====================
function f1 = B_rhs(t,Y,params)
    B=Y(1);
    Q=Y(2);
    X=Y(3);
    r = params(1);
    h=params(2);
    f=params(3);
    g=params(4);
    a=params(5);
    s=params(6);
    k=params(7);
    d=params(8);
    f1 = r*B+g*Q-(h*B*X+f*B);
end

function f2 = Q_rhs(t,Y,params)
    B=Y(1);
    Q=Y(2);
    X=Y(3);
    r=params(1);
    h=params(2);
    f=params(3);
    g=params(4);
    a=params(5);
    s=params(6);
    k=params(7);
    d=params(8);
    f2 = f*B-g*Q;
end


function f3 = X_rhs(t,Y,params)
    B=Y(1);
    Q=Y(2);
    X=Y(3);
    r=params(1);
    h=params(2);
    f=params(3);
    g=params(4);
    a=params(5);
    s=params(6);
    k=params(7);
    d=params(8);
    f3 = a+s*X*(B/(k+B))-d*X;
end
 
function g = rhs(t,Y,params)
    B=Y(1);
    Q=Y(2);
    X=Y(3);
    r=params(1);
    h=params(2);
    f=params(3);
    g=params(4);
    a=params(5);
    s=params(6);
    k=params(7);
    d=params(8);
    g = [B_rhs(t,Y,params);Q_rhs(t,Y,params);X_rhs(t,Y,params)];
end
 



















