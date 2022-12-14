%=====================================================================
% Created on Tue Oct 26 11:34:10 2021
% Noble model using Morris screening with native implementation based on 
% the algorithm in Uncertainty Quantification: Theory, Implementation, 
% and Applications
% by Ralph C. Smith
% Published: 2013
% Pages: XVIII + 382
% Hardcover
% ISBN: 978-1-611973-21-1.
% @author: cogan
%=====================================================================
clear
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');



params=[6,100,120,170,1000,0.1,2,-100,40,-60];
param_labels={'$C_m$', '${\alpha_m}_0$,','${\beta_m}_0$','${\alpha_h}_0$',  '${\beta_h}_0$',  '${\alpha_n}_0$',  '${\beta_n}_0$',   '$E_k$',  '$E_{Na}$',  '$E_{An}$'};
V0=-85;
m0=.5;
h0=.75;
n0=.65;
tspan = linspace(0, 1, 500);

[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)],[V0,m0,h0,n0]); 
%=================================
%Samples
%=================================
p=length(params);
l=20; %number of levels for Morris, larger l means closer approximation, but more samples needed
level=linspace(0,1,l);
delta_star=level(2)-level(1);
Num_samples=200;
d=zeros(Num_samples,p); %Will hold elementary effects
response=zeros(Num_samples,p); % Will hold QoIs for scaling
for i=1:Num_samples
    parameter_matrix=f(params,p,l,level,delta_star);
    for k=1:p
            [t1,y_solution1]=ode45(@(t,Y) rhs(t,Y,parameter_matrix(k,:)), [tspan(1),tspan(end)],[V0,m0,h0,n0]); 
            temp1=abs(real(fft(y_solution1)));
            [M,id1]=max(histcounts(temp1));
            QoI1=temp1(id1);
            [t2,y_solution2]=ode45(@(t,Y) rhs(t,Y,parameter_matrix(k+1,:)), [tspan(1),tspan(end)],[V0,m0,h0,n0]); 
            temp2=abs(real(fft(y_solution2)));
            [M,id2]=max(histcounts(temp2));
            QoI2=temp2(id2);
            response(i,k)=QoI1;
            d(i,k)=abs((QoI2-QoI1)/delta_star); %stores the elementary effects
    end
end
SA=mean(d,1); %unscaled elementary effects
[SA_sort, I] = sort(SA, 'ascend');
bar([1:length(params)],SA_sort/std(std(response(:,I))));
xticks(1:length(params))
xticklabels(param_labels(I))

xlim([0 length(params)+1])
%ylim([0 1])
ylabel('total Sobol index')

%=======================================
%RHS functions
%=======================================
function f1=rhs_V(t,Y,params)
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
    Ek=params(8);
    ENa=params(9);
    EAn=params(10);
    INa=(400000*m^3*h+140)*(V-ENa);
    Ik=1200*n^4*(V-Ek)+(1200*exp((-V-90)/50)+15*exp((V+90)/60))*(V-Ek);
    Ian=75*(V-ENa);
    f1=-(INa+Ik+Ian)/Cm;
end

function f2=rhs_m(t,Y,params)
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
    Ek=params(8);
    ENa=params(9);
    EAn=params(10);
    alpham=alpham_0*(-V-48)/(exp((-V-48)/15)-1);
    betam=betam_0*(V+8)/(exp((V+8)/5)-1);
    f2=alpham*(1-m)-betam*m;
end

function f3=rhs_h(t,Y,params)
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
    Ek=params(8);
    ENa=params(9);
    EAn=params(10);
    alphah=alphah_0*exp((-V-90)/20);
    betah=betah_0/(1+exp((-V-42)/10));
    f3 = alphah*(1-h)-betah*h;
end

    
function f4=rhs_n(t,Y,params)
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
    Ek=params(8);
    ENa=params(9);
    EAn=params(10);
    INa=(400000*m^3*h+140)*(V-ENa);
    alphan=alphan_0*(-V-50)/(exp((-V-50)/10)-1);
    betan=betan_0*exp((-V-90)/80);
    f4 =  alphan*(1-n)-betan*n;
end


function g=rhs(t,Y,params)
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
    Ek=params(8);
    ENa=params(9);
    EAn=params(10);
    INa=(400000*m^3*h+140)*(V-ENa);
    Ik=1200*n^4*(V-Ek)+(1200*exp((-V-90)/50)+15*exp((V+90)/60))*(V-Ek);
    Ian=75*(V-ENa);
    g= [rhs_V(t,Y,params);rhs_m(t,Y,params);rhs_h(t,Y,params);rhs_n(t,Y,params)];
end


function parameters_scaled = f(params,p,l,level,delta_star)
%=======================================================
%Build the Morris method matrix from the 
%algorithm in:
%Uncertainty Quantification Theory, Implementation, 
% and Applications
%By Ralph C. Smith · 2013
%The main idea is to make a matrix that moves randomly
%in the ith parameter direction
%this function gives random steps in each direction to measure 
%the elementary effects
%=======================================================

%====================
%Generate D_star
%====================
rand_zero=randi([0 1], 1, p); %generates random array of 0s and 1s
rand_zero(find(rand_zero==0))=-1; % moves this to -1s and 1s
D_star=diag(rand_zero);
%====================
%Generate P_star
%====================
P_star=eye(p);
P_star(:, randperm(size(P_star, 2)));
%====================
%Generate B
%====================
B=tril(ones(p,p),-1);
B=[B;ones(1,p)];
%====================
%Generate J
%====================
J=ones(p+1,p);
%====================
%Random seed in unit hypercube
%====================
id_rand=randi([1 10],p,1);
q_star=level(id_rand);
%====================
%parameter sample matrix
%====================
B_star=(J*q_star'+delta_star/2*((2*B-J)*D_star+J))*P_star;
%====================
%Scale to get parameters
%====================
l_bounds = .95*params;
u_bounds = 1.05*params;
parameters_scaled=zeros(size(B_star));
    for i=1:p+1
        parameters_scaled(i,:)=(u_bounds-l_bounds).*B_star(i,:)+l_bounds;
    end
end
