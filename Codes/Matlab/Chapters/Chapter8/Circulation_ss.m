

% =============================================================================
% Created on Tue Oct 26 11:34:10 2021
% Circulation model: Algebraic models demonstrating sampling methods
% @author: cogan
% =============================================================================

clear
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
options = optimoptions(@fsolve, 'MaxFunctionEvaluations', 10000, 'Display', 'off');
% =============================================================================
% We can find at least one root
% =============================================================================
params=[5, 5.5, 1.50, 11,2,80, 5,4, .5,.12];
solution1 = fsolve(@(Y)rhs(Y,params), [10,.11,.1],options)
% =============================================================================
% Normal distribution
% =============================================================================
n_power=10;
Num_samples=2^n_power;
total_bacteria=zeros(1,Num_samples);


parameters_normal=lhsnorm(params,diag(.01*params),Num_samples);
% =============================================================================
% Uniform distribution
% =============================================================================
% =============================================================================
% Note that lhsdesign is defined on the unit interval
% =============================================================================
parameters_uniform_unscaled = lhsdesign(Num_samples,length(params));
% =============================================================================
% Scale the samples into the correct parameter scale
% =============================================================================
l_bounds = .95*params;
u_bounds = 1.05*params;
parameters_uniform=zeros(size(parameters_uniform_unscaled));
for i=1:Num_samples
    parameters_uniform(i,:)=(u_bounds-l_bounds).*parameters_uniform_unscaled(i,:)+l_bounds;
end
% =============================================================================
% Sobol distribution: Sobol sequences are low discrepancy sequence
% sobolset is a Matlab native function with many variants -- this is a
% common one
% =============================================================================
p = sobolset(length(params),'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
parameters_sobol_unscaled=net(p,Num_samples);
% =============================================================================
% Scale the samples into the correct parameter scale
% =============================================================================
l_bounds = .95*params;
u_bounds = 1.05*params;
parameters_sobol=zeros(size(parameters_sobol_unscaled));
for i=1:Num_samples
    parameters_sobol(i,:)=(u_bounds-l_bounds).*parameters_sobol_unscaled(i,:)+l_bounds;
end
FP_uniform=zeros(3,Num_samples);
FP_normal=zeros(3,Num_samples);
FP_sobol=zeros(3,Num_samples);


for i=1:Num_samples
    [FP_uniform(:,i),flag1]=fsolve(@(Y)rhs(Y,parameters_uniform(i,:)), [60,2.5,.1],options);
    [FP_normal(:,i),flag2]=fsolve(@(Y)rhs(Y,parameters_normal(i,:)), [60,2.5,.1],options);
    [FP_sobol(:,i),flag3]=fsolve(@(Y)rhs(Y,parameters_sobol(i,:)), [60,2.5,.1],options);
end
figure()
histogram(FP_uniform(1,:),20);
hold on
histogram(FP_normal(1,:),20);
histogram(FP_sobol(1,:),20);
xlabel('$QoI=Q$', fontsize=16)
ylabel('Frequency', fontsize=16)
legend('Uniform', 'Normal', 'Sobol')

figure()
histogram(FP_uniform(2,:),20);
hold on
histogram(FP_normal(2,:),20);
histogram(FP_sobol(2,:),20);
xlabel('$QoI=Q$', fontsize=16)
ylabel('Frequency', fontsize=16)
legend('Uniform', 'Normal', 'Sobol')

figure()
histogram(FP_uniform(3,:),20);
hold on
histogram(FP_normal(3,:),20);
histogram(FP_sobol(3,:),20);
xlabel('$QoI=Q$', fontsize=16)
ylabel('Frequency', fontsize=16)
legend('Uniform', 'Normal', 'Sobol')

% =============================================================================
% Define the RHS of ODEs
% =============================================================================


function  f1 = rhs_Q(Y,params)
    Q=Y(1);
    Pa=Y(2);
    Pv=Y(3);
    Qvessel=params(1);
    R=params(2);
    gamma=params(3);
    Vvessel=params(4);
    V0=params(5);
    F=params(6);
    Vmax=params(7);
    Vmin=params(8);
    Cd=params(9);
    Cs=params(10);
    f1 = (Pa-Pv)/R*(1+gamma*(Pa+Pv)+gamma^2/3*(Pa^2+Pa*Pv+Pv^2))-Qvessel;
end

function  f2 = rhs_Pa(Y,params)
    Q=Y(1);
    Pa=Y(2);
    Pv=Y(3);
    Qvessel=params(1);
    R=params(2);
    gamma=params(3);
    Vvessel=params(4);
    V0=params(5);
    F=params(6);
    Vmax=params(7);
    Vmin=params(8);
    Cd=params(9);
    Cs=params(10);
    f2 = V0*(1+gamma/2*(Pa+Pv)+gamma^2/6*(Pa-Pv)^2)-Vvessel;
end

function  f3 = rhs_Pv(Y,params)
    Q=Y(1);
    Pa=Y(2);
    Pv=Y(3);
    Qvessel=params(1);
    R=params(2);
    gamma=params(3);
    Vvessel=params(4);
    V0=params(5);
    F=params(6);
    Vmax=params(7);
    Vmin=params(8);
    Cd=params(9);
    Cs=params(10);
    f3 = F*(Vmax-Vmin+Cd*Pv-Cs*Pa)-Q;
end

function  g = rhs(Y,params)
    Q=Y(1);
    Pa=Y(2);
    Pv=Y(3);
    Qvessel=params(1);
    R=params(2);
    gamma=params(3);
    Vvessel=params(4);
    V0=params(5);
    F=params(6);
    Vmax=params(7);
    Vmin=params(8);
    Cd=params(9);
    Cs=params(10);
    g = [rhs_Q(Y,params),rhs_Pa(Y,params),rhs_Pv(Y,params)];
end

   