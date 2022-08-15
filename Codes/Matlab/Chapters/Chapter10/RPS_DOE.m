%====================================================================
% Created on Tue Oct 26 11:34:10 2021
% Rock-Paper-scissors model. DOE model SA based on the paper:
%    Sarah C Cotter. A screening design for factorial experiments with
%    interactions. Biometrika, 66(2):317–320, 1979. 
% @author: cogan
%====================================================================


clear
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
params_list={'$\alpha_1$','$\alpha_2$','$\alpha_3$', '$\beta_1$',  '$\beta_2$',  '$\beta_3$'};
params=[.2,.2,.2,2,2,2];
Rr0=.33;
Rp0=.1;
Rs0=.1;
tspan = linspace(0, 200, 4000);
[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)],[Rr0,Rp0,Rs0]); 
figure(1)
plot(t,y_solution(:,1),'b','LineWidth',2)
hold on
plot(t,y_solution(:,2),'k','LineWidth',2)
plot(t,y_solution(:,3),'c','LineWidth',2)

xlabel('Time', fontsize=16)
ylabel('Population', fontsize=16)
legend('$P_r$', '$P_p$','$P_s$', fontsize=16)

figure(2)
 
 params=[.025,.025,.025,4,4,4];
 Rr0=.2;
 Rp0=.1;
 Rs0=.3;
tspan = linspace(0, 100, 4000);
[t1,y_solution1]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)],[Rr0,Rp0,Rs0]); 
plot(t1,y_solution1(:,1),'b','LineWidth',2)
hold on
plot(t1,y_solution1(:,2),'k','LineWidth',2)
plot(t1,y_solution1(:,3),'c','LineWidth',2)
% =============================================================================
% make the design matrix indicating high and low values
% Matlab has the command ff2n which is a 2 level full factorial. Matlab
%defaults to levels between 0 and 1. To match the design in the Python
%code and the Cotter paper, we need to scale this.
% =============================================================================


DOE=1-2*ff2n(length(params)); 

params=[.2,.2,.2,2,2,2];
Rr0=.1;
Rp0=.1;
Rs0=.1;

responser=zeros(1,2^length(params));     %do this for each species
responsep=zeros(1,2^length(params));    
responses=zeros(1,2^length(params));    

for k = 1:2^length(params)
    params1=params+.9*params.*DOE(k,:);
    %To get a regular total time, interpolate onto a regular grid
    [t,y_solution]=ode45(@(t,Y) rhs(t,Y,params1), [tspan(1),tspan(end)],[Rr0,Rp0,Rs0]); 
    y_interp=interp1(t,y_solution,tspan);
    responser(k)=sum(tspan(find(y_interp(:,1)>.8)))/tspan(end);
    responsep(k)=sum(tspan(find(y_interp(:,2)>.8)))/tspan(end);
    responses(k)=sum(tspan(find(y_interp(:,3)>.8)))/tspan(end);

end


Cotter_SAr=zeros(1,length(params));
Cotter_SAp=zeros(1,length(params));
Cotter_SAs=zeros(1,length(params));

for i = 1:length(params)
    Cotter_SAr(length(params)-i+1)=sum(DOE(:,i).*responser')/2^length(params)/2;
    Cotter_SAp(length(params)-i+1)=sum(DOE(:,i).*responsep')/2^length(params)/2;
    Cotter_SAs(length(params)-i+1)=sum(DOE(:,i).*responses')/2^length(params)/2;
end
bar_width=.3
figure()
bar([1:length(params_list)]-1*bar_width,Cotter_SAr, bar_width)
hold on
bar([1:length(params_list)]-0*.5*bar_width,Cotter_SAp, bar_width)
bar([1:length(params_list)]+1*bar_width,Cotter_SAs, bar_width)
legend('$P_r$', '$P_p$', '$P_s$',fontsize=16)
xticks([1:length(params_list)])
xticklabels(params_list)
ylabel('SA',fontsize=16)
xlabel('Parameter',fontsize=16)

%=======================================
% RHS functions
%=======================================

function f1 = rhs_Pr(t,Y,params)
    Pr=Y(1);
    Pp=Y(2);
    Ps=Y(3);
    alpha_1=params(1);
    alpha_2=params(2);
    alpha_3=params(3);
    beta_1=params(4);
    beta_2=params(5);
    beta_3=params(6);
    f1 =   Pr*(1-Pr-alpha_1*Pp-beta_1*Ps);
end

function f2 = rhs_Pp(t,Y,params)
    Pr=Y(1);
    Pp=Y(2);
    Ps=Y(3);
    alpha_1=params(1);
    alpha_2=params(2);
    alpha_3=params(3);
    beta_1=params(4);
    beta_2=params(5);
    beta_3=params(6);
    f2 =  Pp*(1-Pp-alpha_2*Ps-beta_2*Pr);
end

function f3 = rhs_Ps(t,Y,params)
    Pr=Y(1);
    Pp=Y(2);
    Ps=Y(3);
    alpha_1=params(1);
    alpha_2=params(2);
    alpha_3=params(3);
    beta_1=params(4);
    beta_2=params(5);
    beta_3=params(6);
    f3 =  Ps*(1-Ps-alpha_1*Pr-beta_1*Pp);
end

function g = rhs(t,Y,params)
    Pr=Y(1);
    Pp=Y(2);
    Ps=Y(3);
    alpha_1=params(1);
    alpha_2=params(2);
    alpha_3=params(3);
    beta_1=params(4);
    beta_2=params(5);
    beta_3=params(6);
    g = [rhs_Pr(t,Y,params);rhs_Pp(t,Y,params);rhs_Ps(t,Y,params)];
end
