
% We use a generic numerical estimation for the jacobian
% based on finite differences. There are more sophisticated ways, but this 
% is reasonably effective. 
clear
close all
%Parameters from the paper

params=[.1181,.3743,1.131,20.19,3.11*10^(-3),1.636,.5*10^3];
%Find the steady-states-we adjust the starting 
%values to find them.
%Because the scales are quite different we introduce different scalings for
%each variable
%fixed_points(E_min,E_max, delta_E, Tol_E, T_min,T_max, delta_T,Tol_T,params)
fp1=fixed_points(0,1,.001,.001,0,1,.001,.001,params);
%Note that this requires a bit of tuning to get this correct
%also illustrating difficulties in generalizing things
%here fp often returns multiple values -- these are all
%'close' to each other. For now, we will take the mean
fp1=mean(fp1);

fp2=fixed_points(0,2,.001,.001,7,10,.02,.01,params);
fp2=mean(fp2);

fp3=fixed_points(0,1,.001,.001,200,300,.5,.1,params);
fp3=mean(fp3);

fp4=fixed_points(0,1,.001,.001,400,600,.1,.05,params);
fp4=mean(fp4);


%Calculate the Jacobian
J = Myjac(0,fp4,params,.001);
%Estimate the eigenvalues
E_vals=eigs(J);



function f=rhs_T(t,Y,params)
    T=Y(1);
    E=Y(2);
    s=params(1);
    d=params(2);
    p=params(3);
    g=params(4);
    m=params(5);
    r=params(6);
    k=params(7);
    f=s-d*E+p*E*T/(g+T)-m*E*T;
end

function g=rhs_E(t,Y,params)
    T=Y(1);
    E=Y(2);
    s=params(1);
    d=params(2);
    p=params(3);
    g=params(4);
    m=params(5);
    r=params(6);
    k=params(7);
    g =r*T*(1-T/k)-E*T;
end

function [fp]= fixed_points(E_min,E_max, delta_E, Tol_E, T_min,T_max, delta_T,Tol_T,params)
fp=[];
i=0;

for x=T_min:delta_T:T_max
      i=i+1;
      j=0;
    for y=E_min:delta_E:E_max
        j=j+1;
        if (abs(rhs_T(0,[x,y],params))<=Tol_T)&&(abs(rhs_E(0,[x,y],params))<=Tol_E)
            fp=[fp;[y,x]];
        end
    end
end
end


%A very simple way to build the Jacobian
function J = Myjac(t,Y,params,epsilon)
    T=Y(1);
    E=Y(2);
    f_E = (rhs_T(0,[E+epsilon,T],params)-rhs_T(0,[E-epsilon,T],params))/epsilon/2;
    f_T = (rhs_T(0,[E,T+epsilon],params)-rhs_T(0,[E,T-epsilon],params))/epsilon/2;
    g_E = (rhs_E(0,[E+epsilon,T],params)-rhs_E(0,[E-epsilon,T],params))/epsilon/2;
    g_T = (rhs_E(0,[E,T+epsilon],params)-rhs_E(0,[E,T-epsilon],params))/epsilon/2;
    J = [[f_T,f_E];[g_T,g_E]];
end

