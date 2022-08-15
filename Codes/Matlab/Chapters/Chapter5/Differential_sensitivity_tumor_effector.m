
%Differential Sensitivity for the 
%Tumor/Effector system

clear
close all

%Same as usual: define f and g, rhs returns the right hand side
%parameters are from the paper, beware the order though

params=[.1181,.3743,1.131,20.19,3.11*10^(-3),1.636,.5*10^3];
%T=Tumor cells
%E=effector cells 
tspan = linspace(0, 50.000, 200);

% =============================================================================
% There are many different sensitivies we can explore
%There are four different steady-states and 7 parameters
%Here SA_1, SA_2, SA_3, and SA_4 denote the sensitivities near
%different steady-states. Each for loop calculates the solution
%to the augmented system
% =============================================================================
SA_1=zeros(7,2);
for j =1:7
    [t,y_SA1]=ode45(@(t,Y) rhs_diff(t,Y,params,j), [0 100],[.1,.10,0,0]);
    SA_1(j,:)=y_SA1(end,3:4);
end

SA_2=zeros(7,2);
for j =1:7
    [t,y_SA2]=ode45(@(t,Y) rhs_diff(t,Y,params,j), [0 100],[.1,10,0,0]); 
    SA_2(j,:)=y_SA2(end,3:4);
end

SA_3=zeros(7,2);
for j =1:7
    [t,y_SA3]=ode45(@(t,Y) rhs_diff(t,Y,params,j), [0 100],[1,280,0,0]); 
    SA_3(j,:)=y_SA3(end,3:4);
end

SA_4=zeros(7,2);
for j =1:7
    [t,y_SA4]=ode45(@(t,Y) rhs_diff(t,Y,params,j), [0 100],[0.001, 440 ,0,0]); 
    SA_4(j,:)=y_SA4(end,3:4);
end

    
% =============================================================================
% These can be plotted using a bar plot to compare
%how the tumor and effector QoIs compare
% =============================================================================
bar((1:7),[SA_4(:,1),SA_4(:,2)])
xticks(1:7)
xticklabels({'$s$', '$d$', '$p$', '$g$', '$m$', '$r$', '$k$'}, 'Interpreter','latex')
% %  =============================================================================
% %  To investigate how the direct sensitivities change in time for each parameter
% % we can look at the solution to the system
% %  =============================================================================
figure()
j=2
[t,y_time]=ode45(@(t,Y) rhs_diff(t,Y,params,j), [0 100],[1,10 ,0,0]); 
plot(t,y_time(:,3:4),LineWidth=2)
xlabel('Time')
ylabel('Y')
legend('Effector','Tumor')



function f = rhs_E(t,Y,params)
    E=Y(1);
    T=Y(2);
    s=params(1);
    d=params(2);
    p=params(3);
    g=params(4);
    m=params(5);
    r=params(6);
    k=params(7);
    f = s-d*E+p*E*T/(g+T)-m*E*T;
end

function g = rhs_T(t,Y,params)
    E=Y(1);
    T=Y(2);
    s=params(1);
    d=params(2);
    p=params(3);
    g=params(4);
    m=params(5);
    r=params(6);
    k=params(7);
    g =r*T*(1-T/k)-E*T;
end
function h = rhs(t,Y,params)
    E=Y(1);
    T=Y(2);
    s=params(1);
    d=params(2);
    p=params(3);
    g=params(4);
    m=params(5);
    r=params(6);
    k=params(7);
    h = [f(t,Y,params);g(t,Y,params)];
end

% j will denote the parameter we are looking at
function A = rhs_diff(t,Y,params,j)
%Augment with s1, s2 differential for T and E
    E=Y(1);
    T=Y(2);
    s1=Y(3);
    s2=Y(4);
    s=params(1);
    d=params(2);
    p=params(3);
    g=params(4);
    m=params(5);
    r=params(6);
    k=params(7);
%increment the parameter
    dp=.01*params(j);
    dE=.001; %for e component of Jac
    dT=.001; %for t component of Jac
    params_p= params;
    params_m=params;
    params_p(j)=params(j)+dp;
    params_m(j)=params(j)-dp;
    f1 = (rhs_E(0,[E,T],params_p)-rhs_E(0,[E,T],params_m))/(2*dp);
    g1 = (rhs_T(0,[E,T],params_p)-rhs_T(0,[E,T],params_m))/(2*dp);
    f_e = (rhs_E(0,[E+dE,T],params)-rhs_E(0,[E-dE,T],params))/(2*dE);
    f_t = (rhs_E(0,[E,T+dT],params)-rhs_E(0,[E,T-dT],params))/(2*dT);
    g_e = (rhs_T(0,[E+dE,T],params)-rhs_T(0,[E-dE,T],params))/(2*dE);
    g_t = (rhs_T(0,[E,T+dT],params)-rhs_T(0,[E,T-dT],params))/(2*dT);
    %[f_e,f_t;g_e,g_t] is the Jac
    A = [rhs_E(0,[E,T],params);rhs_T(0,[E,T],params);f1+f_e*s1+f_t*s2;g1+g_e*s1+g_t*s2];
end
