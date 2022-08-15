%================================================================
% Created on Tue Oct 26 11:34:10 2021
%Hodgkin-Huxley model: Sobol Sensitivity as implemented in SAlib
% @author: cogan
% version 1 uses: http://codes.arizona.edu/toolbox/help/html/installation.html
% http://codes.arizona.edu/toolbox/help/html/CODES.html
% Download the CODES package
% Need to run the following commands:
% addpath('CODES/')  
% CODES.install
% We also include native code written by: Alen Alexanderian, Pierre Gremaud, and Ralph Smith
% https://aalexan3.math.ncsu.edu/rtg/rtg_slides.pdf
%================================================================


clear
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');



params_list=[ '$C_m$', '${\alpha_m}_0$', '${\beta_m}_0$', '${\alpha_h}_0$', '${\beta_h}_0$', '${\alpha_n}_0$', '${\beta_n}_0$', '$g_K$', '$g_{Na}$', '$g_l$', '$V_k$', '$V_{Na}$', '$V_l$', '$I_{app}$'];

%===============================
%Compare different 
%behaviors
%===============================

params1=[2.000, .1, 4,.07,1,0.01,.125,36, 120, .3,-12, 115,10.6,0*35];
V0=1;
m0=.053;
h0=.595;
n0=.317;
tspan = linspace(0, 20, 500);
[t1,y_solution1]=ode45(@(t,Y) rhs(t,Y,params1), [tspan(1),tspan(end)],[V0,m0,h0,n0]); 

params2=[2.000, .1, 4,.07,1,0.01,.125,36, 120, .3,-12, 115,10.6,0*35];
V0=6;
m0=.053;
h0=.595;
n0=.317;
tspan = linspace(0, 20, 500);
[t2,y_solution2]=ode45(@(t,Y) rhs(t,Y,params2), [tspan(1),tspan(end)],[V0,m0,h0,n0]); 

params3=[2.000, .1, 4,.07,1,0.01,.125,36, 120, .3,-12, 115,10.6,0*35];
V0=20;
m0=.053;
h0=.595;
n0=.317;
[t3,y_solution3]=ode45(@(t,Y) rhs(t,Y,params3), [tspan(1),tspan(end)],[V0,m0,h0,n0]); 

plot(t1,y_solution1(:,1),'b','LineWidth',2)
hold on
plot(t2,y_solution2(:,1),'m','LineWidth',2)
plot(t3,y_solution3(:,1),'k','LineWidth',2)
legend('$V_0=1$','$V_0=6$','$V_0=20$')
xlabel('Time', fontsize=16)
ylabel('Voltage', fontsize=16)
% =============================================================================
% Sobol' sensitivity. We are using the CODES package since native
% Sobol' is complex and often far to slow.
% =============================================================================

% =============================================================================
% First CODES requires simple input/output structure. We define a new
% function SA =QoI(params). The first example gives the maximum value 
% of the voltage for a given set of parameters, initial conditions and
% time-span
% =============================================================================

Test=QoI(params1) 


dim=14;
n=1e5;
%%Sobol_out=CODES.sensitivity.sobol(@(params1) QoI(params1),dim,n,'bar_plot',true,'conv_seq',linspace(100,n,20));
%%Sobol_out.S1 %gives a table with indices


% =============================================================================
% Sobol' sensitivity using codes developed by Ralph Smith
% =============================================================================

%
% Computes the Sobol' indices for the scalar QoI
%
ndim = length(params1);
Ns = 1e4;

for i = 1 : ndim
   fprintf('param %i analysis \n', i);
   [Xi F1, F2, mu, D, Stot(i)] = get_sobol_scalar(@QoI, ndim, Ns, i, params1);
end

close all;
figure(1);
param_labels = {'$C_m$', '${\alpha_m}_0$', '${\beta_m}_0$', '${\alpha_h}_0$', '${\beta_h}_0$', '${\alpha_n}_0$', '${\beta_n}_0$', '$g_K$', '$g_{Na}$', '$g_l$', '$V_k$', '$V_{Na}$', '$V_l$', '$I_{app}$'}; 
bar(Stot);
set(gca,'fontsize', 20, 'xticklabels', param_labels); 
xlim([0 length(params1)])
ylim([0 1])
ylabel('total Sobol index')

%================================================================
% Define the QoI function Sobol'
%================================================================
function SA = QoI(params)
    V0=6;
    m0=.053;
    h0=.595;
    n0=.317;
    tspan = linspace(0, 20, 500);
 %====================
     [t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)],[V0,m0,h0,n0]); 
     SA = max(y_solution(:,1));
end


%================================================================
%Beginning of Scalar Sobol' function
%================================================================

function [Xi F1, F2, mu sigma2 Sy_tot Sy] = get_sobol_scalar(F, ndim, Ns, yidx, params1)
%
%  F    : function handle to simulation routine
%         the goal is to assess sensitivity of F(x_1, ..., x_ndim) to parameters, x_1, ..., x_ndim
%
%  ndim : stochastic dimension 
%
%  Ns   : the number of samples used to estimate the sensitivity indices
%
%  yidx : the index set for which we want to compute the index,
%
%         example yidx = [1 2] gives S_12, yidx = [2] gives S_2, etc.
%         the following is formulated following the paper of Sobol ...
% 
% reference: 
%   Sobol, I.M., Global sensitivity indices for nonlinear mathematical models and their Monte Carlo estimates, 
%   Mathematics and Computers in Simulation, 2001.

m = length(yidx);
p = ndim - m;

% get the basic samples needed
eta        = 2*rand(Ns, m)-1;
etaprime   = 2*rand(Ns, m)-1;
zeta       = 2*rand(Ns, p)-1;
zetaprime  = 2*rand(Ns, p)-1;

% index sets
I = [1 : ndim];     
idx = setdiff(I, yidx);

% set up the samples
Xi1 = zeros(Ns, ndim);
Xi2 = zeros(Ns, ndim);
Xi3 = zeros(Ns, ndim);

Xi1(:, yidx) = eta;
Xi1(:, idx)  = zeta;

Xi2(:, yidx)   = eta;
Xi2(:, idx) = zetaprime;

Xi3(:, yidx)   = etaprime;
Xi3(:, idx) = zeta;


% get the needed function evaluations 
F1 = zeros(Ns, 1);
F2 = zeros(Ns, 1);
F3 = zeros(Ns, 1);
% the sampling loop
for k = 1 : Ns
   ratio = (k/Ns)*100;
       if mod(ratio, 5)<1e-13
           if ratio == 5
              fprintf('[');
           end
       fprintf('%g%% ', ratio);
           if ratio == 100
              fprintf(']\n');
           end
        end

%Note that this code is set-up so that Xi1 provides an index to parameters
%we need to scale these to true parameter ranges. 
% In the original code, this is done in the ODEsolver
% To maintain the structure here I am sacrificing some 
%efficiency 
    meanpar=params1;
    a = meanpar - 0.2 * meanpar;
    b = meanpar + 0.2 * meanpar;
    paramsi1  = 0.5 * (a + b) + 0.5 * (b - a) .* Xi1(k, :);
    paramsi2  = 0.5 * (a + b) + 0.5 * (b - a) .* Xi2(k, :);
    paramsi3  = 0.5 * (a + b) + 0.5 * (b - a) .* Xi3(k, :);

   F1(k,:) = F(paramsi1);
   F2(k,:) = F(paramsi2);
   F3(k,:) = F(paramsi3);
end

% Compute the estimators
phi  = F1;
phi2 = F1.^2;
psi  = phi .* F2;
chi  = 0.5 * (phi - F3).^2;

% compute the results
mu      = mean(phi);
mu2     = mean(phi2);
sigma2  = mu2 - mu.^2;
mupsi   = mean(psi);
Sy      = (mupsi - mu.^2) ./ sigma2;
Sy_tot  = mean(chi) ./ sigma2;

Xi = Xi1;

mu     = mu(:);
sigma2 = sigma2(:);
Sy     = Sy(:);
Sy_tot = Sy_tot(:);
end

%================================================================
%End of Scalar Sobol' function
%================================================================
%================================================================
% Define the RHS functions
%================================================================


function f1 = rhs_V(t,Y,params)
    V=Y(1);
    m=Y(2);
    h=Y(3);
    n=Y(4);
    Cm=params(1);
    alpham_0=params(2);
    betam_0=params(3);
    alphah_0=params(4);
    betah_0=params(5);
    alphan_0=params(6);
    betan_0=params(7);
    gK=params(8);
    gNa=params(9);
    gl=params(10);
    Vk=params(11);
    VNa=params(12);
    Vl=params(13);
    Iapp=params(14);
    INa=gNa*m^3*h*(V-VNa);
    Ik=gK*n^4*(V-Vk);
    Il=gl*(V-Vl);
    f1 = -(INa+Ik+Il)/Cm+Iapp;
end

    
function f2 = rhs_m(t,Y,params)
    V=Y(1);
    m=Y(2);
    h=Y(3);
    n=Y(4);
    Cm=params(1);
    alpham_0=params(2);
    betam_0=params(3);
    alphah_0=params(4);
    betah_0=params(5);
    alphan_0=params(6);
    betan_0=params(7);
    gK=params(8);
    gNa=params(9);
    gl=params(10);
    Vk=params(11);
    VNa=params(12);
    Vl=params(13);
    alpham=alpham_0*(25-V)/(exp((25-V)/10)-1);
    betam=betam_0*exp(-V/18);
    f2 = alpham*(1-m)-betam*m;
end

function f3 = rhs_h(t,Y,params)
    V=Y(1);
    m=Y(2);
    h=Y(3);
    n=Y(4);
    Cm=params(1);
    alpham_0=params(2);
    betam_0=params(3);
    alphah_0=params(4);
    betah_0=params(5);
    alphan_0=params(6);
    betan_0=params(7);
    gK=params(8);
    gNa=params(9);
    gl=params(10);
    Vk=params(11);
    VNa=params(12);
    Vl=params(13);
     alphah=alphah_0*exp((-V)/20);
    betah=betah_0/(1+exp((-V+30)/10));
    f3 = alphah*(1-h)-betah*h;
end


function f4 = rhs_n(t,Y,params)
    V=Y(1);
    m=Y(2);
    h=Y(3);
    n=Y(4);
    Cm=params(1);
    alpham_0=params(2);
    betam_0=params(3);
    alphah_0=params(4);
    betah_0=params(5);
    alphan_0=params(6);
    betan_0=params(7);
    gK=params(8);
    gNa=params(9);
    gl=params(10);
    Vk=params(11);
    VNa=params(12);
    Vl=params(13);
     alphan=alphan_0*(-V+10)/(exp((-V+10)/10)-1);
    betan=betan_0*exp((-V)/80);
    f4 = alphan*(1-n)-betan*n;
end


function g = rhs(t,Y,params)
    V=Y(1);
    m=Y(2);
    h=Y(3);
    n=Y(4);
    Cm=params(1);
    alpham_0=params(2);
    betam_0=params(3);
    alphah_0=params(4);
    betah_0=params(5);
    alphan_0=params(6);
    betan_0=params(7);
    gK=params(8);
    gNa=params(9);
    gl=params(10);
    Vk=params(11);
    VNa=params(12);
    Vl=params(13);
    g = [rhs_V(t,Y,params);rhs_m(t,Y,params);rhs_h(t,Y,params);rhs_n(t,Y,params)];
end
