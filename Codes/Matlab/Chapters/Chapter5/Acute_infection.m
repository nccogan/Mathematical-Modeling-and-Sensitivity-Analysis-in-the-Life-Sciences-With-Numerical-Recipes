%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Sun Oct 28 22:28:54 2018
% 
% @author: cogan
% Feature SA using post-processing
% Acute infection model
% T: Target cells
% I: Infected cells
% V: Virus cells
% A: antibodies (signals)
% params=[beta, delta, p, c, c_a, k, S_a, d]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

% =============================================================================
% Direct numerical simulation
% =============================================================================

params=[ .1181,.0743,1.131,20.19,3.11*10^(-3),1.636,.5*10^3,1];
tstart=0;
tfinal=150;
N=500;
tspan = linspace(0, 150, 500);
dt=(tfinal-tstart)/(N-1);
y0 = [100.35,1,0,0];
[t,y_simulation]=ode45(@(t,Y) rhs(t,Y,params), [0 100],y0); 
 
figure(1)
plot(t,y_simulation(:,1:2),LineWidth=2)
xlabel('Time', fontsize=16)
ylabel('Y', fontsize=16)
legend('Target','Infected', fontsize=16)

figure(2)
plot(t,y_simulation(:,2),LineWidth=2)
legend(['Virus'], fontsize=16)
xlabel('Time', fontsize=16)
ylabel('Viral Concentration', fontsize=16)


figure(3)
plot(t,y_simulation(:,3),LineWidth=2)
xlabel('Time', fontsize=16)
ylabel('Antibody Concentration', fontsize=16)


%Post-processing to estimate derivative using the gradient operator
figure(4)
y_simulation_interp=interp1(t,y_simulation(:,3),tspan); %interpolate onto regular grid
%Note that here we are not using the time vector defined by ode45
%try it and see what happens!
gradient_v=gradient(y_simulation_interp,dt);
plot(tspan, gradient_v,LineWidth=2)
xlabel('Time', fontsize=16)
ylabel('$\frac{dV}{dt}$', 'Interpreter','latex',fontsize=16)


% =============================================================================
% Feature Sensitivity using post-processing
% =============================================================================
tspan = linspace(tstart, tfinal, N);
y0_diff=[100.35,1,0,0,0,0,0,0];
[t,y_simulation]=ode45(@(t,Y) rhs_diff(t,Y,params,6), [0 100],y0_diff); 

figure(5)
plot(t,y_simulation(:,6),LineWidth=2)
xlabel('Time', fontsize=16)
ylabel('Sens. of $V$ to $S_A$',  'Interpreter','latex',fontsize=16)
% 
% 
% gradient_v_feature=np.gradient(y_SA_feature.y[6],(tspan[2]-tspan[1]));
% 
% pyplot.figure(4)
% pyplot.plot(tspan, gradient_v_feature)
% pyplot.xlabel('Time', fontsize=16)
% pyplot.ylabel('Sens. of ' r'$\frac{\partial V}{\partial t}$ to $S_A$', fontsize=16)



%=================
%Functions
%=================

% =============================================================================
% Define the RHS of ODEs
% This section considers the direct
% simulation.  
% =============================================================================
function f1 = rhs_T(t,Y,params)
    T=Y(1);
    I=Y(2);
    V=Y(3);
    A=Y(4);
    f1 = -params(1)*T*V;
end

function f2 = rhs_I(t,Y,params)
    T=Y(1);
    I=Y(2);
    V=Y(3);
    A=Y(4);
    f2 = params(1)*T*V-params(2)*I;
end

function f3 = rhs_V(t,Y,params)
    T=Y(1);
    I=Y(2);
    V=Y(3);
    A=Y(4);
    f3 = params(3)*I-params(4)*V-params(5)*V*A;
end

function f4 = rhs_A(t,Y,params)
    T=Y(1);
    I=Y(2);
    V=Y(3);
    A=Y(4);
    f4  = params(6)*V*A+params(7)-params(8)*A;
end

function g = rhs(t,Y,params)
    T=Y(1);
    I=Y(2);
    V=Y(3);
    A=Y(4);
    g = [rhs_T(t,Y,params);rhs_I(t,Y,params);rhs_V(t,Y,params);rhs_A(t,Y,params)];
end

% =============================================================================
% Differential SA
% For this we have augmented the 4 odes with the derivatives with respect to the
%  parameters. We have to indicate whch parameter -- here this is j
% Note that we could be smarter about generalizing our Jacobians rather than typing them by hand 
% =============================================================================


function g_diff =  rhs_diff(t,Y,params,j) % j will denote the parameter we are looking at
    T=Y(1);
    I=Y(2);
    V=Y(3);
    A=Y(4);
    s1=Y(5);
    s2=Y(6);
    s3=Y(7);
    s4=Y(8);
    dp=.1*params(j); %increment the parameter
    dy=.001; %for component of Jac
    params_p= params;
    params_m=params;
    params_p(j)=params(j)+dp;
    params_m(j)=params(j)-dp;
%f_p entries
    fp_1 = (rhs_T(t,[T,I,V,A],params_p)-rhs_T(t,[T,I,V,A],params_m))/(2*dp);
    fp_2 = (rhs_I(t,[T,I,V,A],params_p)-rhs_I(t,[T,I,V,A],params_m))/(2*dp);
    fp_3 = (rhs_V(t,[T,I,V,A],params_p)-rhs_V(t,[T,I,V,A],params_m))/(2*dp);
    fp_4 = (rhs_A(t,[T,I,V,A],params_p)-rhs_A(t,[T,I,V,A],params_m))/(2*dp);
%We could loop through to get the Jacobian but we are doing it explicitly
%First row
    f1_1 = (rhs_T(t,[T+dy,I,V,A],params)-rhs_T(t,[T-dy,I,V,A],params))/(2*dy);
    f1_2 = (rhs_T(t,[T,I+dy,V,A],params)-rhs_T(t,[T,I-dy,V,A],params))/(2*dy);
    f1_3 = (rhs_T(t,[T,I,V+dy,A],params)-rhs_T(t,[T,I,V-dy,A],params))/(2*dy);
    f1_4 = (rhs_T(t,[T,I,V,A+dy],params)-rhs_T(t,[T,I,V,A-dy],params))/(2*dy);
%Second row
    f2_1 = (rhs_I(t,[T+dy,I,V,A],params)-rhs_I(t,[T-dy,I,V,A],params))/(2*dy);
    f2_2 = (rhs_I(t,[T,I+dy,V,A],params)-rhs_I(t,[T,I-dy,V,A],params))/(2*dy);
    f2_3 = (rhs_I(t,[T,I,V+dy,A],params)-rhs_I(t,[T,I,V-dy,A],params))/(2*dy);
    f2_4 = (rhs_I(t,[T,I,V,A+dy],params)-rhs_I(t,[T,I,V,A-dy],params))/(2*dy);
%Third row
    f3_1 = (rhs_V(t,[T+dy,I,V,A],params)-rhs_V(t,[T-dy,I,V,A],params))/(2*dy);
    f3_2 = (rhs_V(t,[T,I+dy,V,A],params)-rhs_V(t,[T,I-dy,V,A],params))/(2*dy);
    f3_3 = (rhs_V(t,[T,I,V+dy,A],params)-rhs_V(t,[T,I,V-dy,A],params))/(2*dy);
    f3_4 = (rhs_V(t,[T,I,V,A+dy],params)-rhs_V(t,[T,I,V,A-dy],params))/(2*dy);
%Fourth row
    f4_1 = (rhs_A(t,[T+dy,I,V,A],params)-rhs_A(t,[T-dy,I,V,A],params))/(2*dy);
    f4_2 = (rhs_A(t,[T,I+dy,V,A],params)-rhs_A(t,[T,I-dy,V,A],params))/(2*dy);
    f4_3 = (rhs_A(t,[T,I,V+dy,A],params)-rhs_A(t,[T,I,V-dy,A],params))/(2*dy);
    f4_4 = (rhs_A(t,[T,I,V,A+dy],params)-rhs_A(t,[T,I,V,A-dy],params))/(2*dy);

g_diff =  [rhs_T(t,[T,I,V,A],params);rhs_I(t,[T,I,V,A],params);rhs_V(t,[T,I,V,A],params);rhs_A(t,[T,I,V,A],params);
            fp_1+f1_1*s1+f1_2*s2+f1_3*s3+f1_4*s4;
            fp_2+f2_1*s1+f2_2*s2+f2_3*s3+f2_4*s4;
            fp_3+f3_1*s1+f3_2*s2+f3_3*s3+f3_4*s4;
            fp_4+f4_1*s1+f4_2*s2+f4_3*s3+f4_4*s4];            
end


