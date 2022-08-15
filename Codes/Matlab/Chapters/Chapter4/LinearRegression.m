%Finding the slopes of the 
%regression curves
clear
close all


%Parameter Descriptions
%params=[r1,kappa1,alpha12,r2,kappa2,alpha21]
%Example parameters for different casts of the 
%competitive exclusion
%case 1 params=[1,1,2,1,2,1]
%case 2 params=[1,2,1,1,1,2]
%case 3
params=[1,2,1,1,3,2];
%case 4 params=[1,3,2,1,2,1]



%Sample parameters
Num_samples=500;
tstart=0;
tstop=26;
%Initial conditions for N1 and N2
N1_0=.1;
N2_0=.1;
%Bound the parameter sample space
params_max=1.1*[1,2,1,1,3,2];
params_min=.9*[1,2,1,1,3,2];
%Initialize QoI
p_vals=zeros(Num_samples,1);
%Initialize different QoIs
QoI1=zeros(Num_samples,1);
QoI2=zeros(Num_samples,1);


%This loop runs through all parameters. 
%The outside loop (k)
%fixes the parameter of interest.
%The inside loop (s) 
%goes through each sample

for k = 1:6
    params=(params_max+params_min)/2;
    for s = 1: Num_samples+1
        params(k)= (params_max(k)-params_min(k)).*rand(1) + params_min(k);
        tspan = linspace(tstart, tstop, 200);
        Y0 = [N1_0, N2_0];
        [t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tstart tstop], Y0); 
        ys=y_solution;
        QoI1(s)=ys(end,1)/(ys(end,1)+ys(end,2));
        QoI2(s)=ys(end,2)/(ys(end,2)+ys(end,2));
        p_vals(s)=params(k);
   
        coeffs=polyfit(p_vals,QoI1,  1);%Polyfit interpolates the particular QoI vs parameter data 
                                        %generated with a polynomial. The degree of the polynomial is an option an the return is the 
                                        %coefficent. 
    end
     scatter(p_vals, QoI1);
    hold on
    ylabel('QoI1', fontsize = 16)
    %xlabel('%s'%pv[k], fontsize = 16)                                   
     hold on
     
     plot(p_vals,coeffs(2)+coeffs(1)*p_vals);
     figure

 
end

%Functions
function f=rhs_N1(t,Y,params)
    y1=Y(1);
    y2=Y(2);
    r1 = params(1);
    kappa1= params(2);
    alpha12= params(3);
    r2= params(4);
    kappa2= params(5);
    alpha21 = params(6);
    f =  r1*y1*(kappa1-y1-alpha12*y2)/kappa1;
end

function g=rhs_N2(t,Y,params)
    y1=Y(1);
    y2=Y(2);
    r1 = params(1);
    kappa1= params(2);
    alpha12= params(3);
    r2= params(4);
    kappa2= params(5);
    alpha21 = params(6);
    g =  r2*y2*(kappa2-y2-alpha21*y1)/kappa2;
end

function r=rhs(t,Y,params)
    y1=Y(1);
    y2=Y(2);
    r1 = params(1);
    kappa1= params(2);
    alpha12= params(3);
    r2= params(4);
    kappa2= params(5);
    alpha21 = params(6);
    r =  [rhs_N1(t,Y,params); rhs_N2(t,Y,params)];
end

