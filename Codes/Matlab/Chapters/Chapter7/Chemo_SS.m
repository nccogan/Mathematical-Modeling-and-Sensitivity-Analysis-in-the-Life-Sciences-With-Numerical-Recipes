% =============================================================================
% Created on Tue Oct 26 11:34:10 2021
% Chemostat model : Steady-state
% @author: cogan
% =============================================================================

clear;
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

N_0=1;
B_0=.05;
params=[1, .05, .25, .5, 1];
solution1=fsolve(@(Y) rhs(Y,params), [.1, .1])
 J=Myjac([solution1(1),solution1(2)],params,.001)

% =============================================================================
% Determine the eigenvalues
% =============================================================================
w=eigs(J)
[t,y_solution]=ode45(@(t,Y) rhs(Y,params), [0 200],[1.10, .1]); 
plot(t,y_solution(:,1),'c','LineWidth',2)
hold on
plot(t,y_solution(:,2),'k','LineWidth',2)
xlabel('Time', fontsize=16)
ylabel('Concentration', fontsize=16)
legend('Nutrient', 'Bacteria', fontsize=16)


% =============================================================================
% Define the RHS of ODEs for nutrient (N) and bacteria (B)
% =============================================================================

function f1=rhs_N(Y,params)
    N=Y(1);
    B=Y(2);
    N0=params(1);
    F=params(2);
    Yield=params(3);
    mu=params(4);
    K_n=params(5);
    f1 = N0*F-1/Yield*mu*N/(K_n+N)*B-F*N;
end

function f2=rhs_B(Y,params)
    N=Y(1);
    B=Y(2);
    N0=params(1);
    F=params(2);
    Yield=params(3);
    mu=params(4);
    K_n=params(5);
    f2 = mu*N/(K_n+N)*B-F*B;
end

function g=rhs(Y,params)
    N=Y(1);
    B=Y(2);
    N0=params(1);
    F=params(2);
    Yield=params(3);
    mu=params(4);
    K_n=params(5);
    g=[rhs_N(Y,params);rhs_B(Y,params)];
end

% =============================================================================
% A discrete estimate of the Jacobian using centered difference. The user should
%consider this with care.
% =============================================================================
function J =  Myjac(Y,params,epsilon)
    N=Y(1);
    B=Y(2);
    f_N = (rhs_N([N+epsilon,B],params)-rhs_N([N-epsilon,B],params))/epsilon/2;
    f_B = (rhs_N([N,B+epsilon],params)-rhs_N([N,B-epsilon],params))/epsilon/2;
    g_N = (rhs_B([N+epsilon,B],params)-rhs_N([N-epsilon,B],params))/epsilon/2;
    g_B = (rhs_B([N,B+epsilon],params)-rhs_N([N,B-epsilon],params))/epsilon/2;
    J = [[f_N,f_B];[g_N,g_B]];
end


