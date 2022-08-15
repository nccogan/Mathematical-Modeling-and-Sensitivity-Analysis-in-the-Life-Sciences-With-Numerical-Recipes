
% =============================================================================
% Created on Tue Oct 26 11:34:10 2021
%Chemostat model and Pearsons Moment Correlation
%@author: cogan
% =============================================================================


clear
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

N_0=0;
B_0=.05;
% =============================================================================
% Coexistence example
% =============================================================================
params=[1, .05, .25, .5, 1];
tspan = linspace(0, 100, 500);
[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)],[N_0, B_0]); 
figure(1)
plot(t,y_solution(:,1),'c','LineWidth',2)
hold on
plot(t,y_solution(:,2),'k','LineWidth',2)
xlabel('Time', fontsize=16)
ylabel('Concentration', fontsize=16)
legend('Nutrient', 'Bacteria', fontsize=16)

% =============================================================================
% Washout example
% =============================================================================
figure(2)
params=[1, .3, .25, .5, 1];
tspan = linspace(0, 100, 500);
[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)],[N_0, B_0]); 
plot(t,y_solution(:,1),'c','LineWidth',2)
hold on
plot(t,y_solution(:,2),'k','LineWidth',2)
xlabel('Time', fontsize=16)
ylabel('Concentration', fontsize=16)
legend('Nutrient', 'Bacteria', fontsize=16)


% =============================================================================
% Pearsons correlation coefficient
% =============================================================================
Num_samples=50;
total_bacteria=zeros(1,Num_samples);

% =============================================================================
% Note that lhsdesign is defined on the unit interval
% =============================================================================
parameters = lhsdesign(Num_samples,length(params));
% =============================================================================
% Scale the samples into the correct parameter scale
% =============================================================================
l_bounds = .95*params;
u_bounds = 1.05*params;
parameters_scaled=zeros(size(parameters));
for i=1:Num_samples
    parameters_scaled(i,:)=(u_bounds-l_bounds).*parameters(i,:)+l_bounds;
end
% =============================================================================
% Evaluate the model for each parameter set
% =============================================================================
tspan = linspace(0, 100, 500);
for i =1:Num_samples
    params_i=parameters_scaled(i,:);
    [t,y_solution]=ode45(@(t,Y) rhs(t,Y,params_i), [tspan(1),tspan(end)],[N_0, B_0]); 
% =============================================================================
% trapz is a relatively standard  implementation of the 
% trapezoidal rule for integration. QoI is total bacterial count
% =============================================================================
    total_bacteria(i)=trapz(t,y_solution(:,2));
end

cc=zeros(1,5);
for j = 1:length(params)
    CC=corrcoef(parameters_scaled(:,j), total_bacteria);
    cc(j)=CC(1,end);
end

% =============================================================================
% Increase the number of samples
% =============================================================================
Num_samples=500;
total_bacteria=zeros(1,Num_samples);

% =============================================================================
% Note that lhsdesign is defined on the unit interval
% =============================================================================
parameters = lhsdesign(Num_samples,length(params));
% =============================================================================
% Scale the samples into the correct parameter scale
% =============================================================================
l_bounds = .95*params;
u_bounds = 1.05*params;
parameters_scaled=zeros(size(parameters));
for i=1:Num_samples
    parameters_scaled(i,:)=(u_bounds-l_bounds).*parameters(i,:)+l_bounds;
end
% =============================================================================
% Evaluate the model for each parameter set
% =============================================================================
tspan = linspace(0, 100, 500);
for i =1:Num_samples
    params_i=parameters_scaled(i,:);
    [t,y_solution]=ode45(@(t,Y) rhs(t,Y,params_i), [tspan(1),tspan(end)],[N_0, B_0]); 
% =============================================================================
% trapz is a relatively standard  implementation of the 
% trapezoidal rule for integration. QoI is total bacterial count
% =============================================================================
    total_bacteria(i)=trapz(t,y_solution(:,2));
end

cc_500=zeros(1,5);
for j = 1:length(params)
    CC=corrcoef(parameters_scaled(:,j), total_bacteria);
    cc_500(j)=CC(1,end);
end

% =============================================================================
% Even more samples
% =============================================================================
Num_samples=1500;
total_bacteria=zeros(1,Num_samples);

% =============================================================================
% Note that lhsdesign is defined on the unit interval
% =============================================================================
parameters = lhsdesign(Num_samples,length(params));
% =============================================================================
% Scale the samples into the correct parameter scale
% =============================================================================
l_bounds = .95*params;
u_bounds = 1.05*params;
parameters_scaled=zeros(size(parameters));
for i=1:Num_samples
    parameters_scaled(i,:)=(u_bounds-l_bounds).*parameters(i,:)+l_bounds;
end
% =============================================================================
% Evaluate the model for each parameter set
% =============================================================================
tspan = linspace(0, 100, 500);
for i =1:Num_samples
    params_i=parameters_scaled(i,:);
    [t,y_solution]=ode45(@(t,Y) rhs(t,Y,params_i), [tspan(1),tspan(end)],[N_0, B_0]); 
% =============================================================================
% trapz is a relatively standard  implementation of the 
% trapezoidal rule for integration. QoI is total bacterial count
% =============================================================================
    total_bacteria(i)=trapz(t,y_solution(:,2));
end

cc_1500=zeros(1,5);
for j = 1:length(params)
    CC=corrcoef(parameters_scaled(:,j), total_bacteria);
    cc_1500(j)=CC(1,end);
end

figure()
bar_width=.35;
p_names=['$N_0$', '$F$', '$Y$', '$\mu$', '$K_n$' ];
hold on
bar([1:length(params)]-1*bar_width,cc,bar_width ,'b')
bar([1:length(params)]-0*bar_width,cc_500, bar_width,'c')
bar([1:length(params)]+1*bar_width,cc_1500, bar_width,'k')
xticks(1:length(params))
xticklabels({'$N_0$', '$F$', '$Y$', '$\mu$', '$K_n$'})
legend('50 Samples', '500 Samples','1500 Samples')


% =============================================================================
% Don't forget to check linearity
% =============================================================================
figure()

scatter(parameters_scaled(:,1),total_bacteria)
xlabel('$N_0$', 'Interpreter','Latex',fontsize=16)
ylabel('QoI'   , fontsize=16)

figure()
scatter(parameters_scaled(:,2),total_bacteria)
xlabel('$F$','Interpreter','Latex', fontsize=16)
ylabel('QoI'   , fontsize=16)

figure()


scatter(parameters_scaled(:,3),total_bacteria)
xlabel('$N_0$','Interpreter','Latex', fontsize=16)
ylabel('QoI'   , fontsize=16)

figure()

scatter(parameters_scaled(:,4),total_bacteria)
xlabel('$\mu$','Interpreter','Latex', fontsize=16)
ylabel('QoI'   , fontsize=16)

figure()
scatter(parameters_scaled(:,5),total_bacteria)
xlabel('$K_n$','Interpreter','Latex', fontsize=16)
ylabel('QoI'   , fontsize=16)



% =============================================================================
% Define the RHS of ODEs
% =============================================================================

function f1 = rhs_N(t,Y,params)
    N=Y(1);
    B=Y(2);
    N0=params(1);
    F=params(2);
    Yield=params(3);
    mu=params(4);
    K_n = params(5);
    f1 = N0*F-1/Yield*mu*N/(K_n+N)*B-F*N;
end

function f2 = rhs_B(t,Y,params)
    N=Y(1);
    B=Y(2);
    N0=params(1);
    F=params(2);
    Yield=params(3);
    mu=params(4);
    K_n = params(5);
    f2=mu*N/(K_n+N)*B-F*B;
end

function g= rhs(t,Y,params)
    N=Y(1);
    B=Y(2);
    N0=params(1);
    F=params(2);
    Yield=params(3);
    mu=params(4);
    K_n = params(5);
    g=[rhs_N(t,Y,params);rhs_B(t,Y,params)];
end


