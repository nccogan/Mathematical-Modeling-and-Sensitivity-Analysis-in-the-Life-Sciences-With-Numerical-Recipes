%basic Scalar ODE 
clear
close all

%define parameters
params=[.1,100];
%Define the time interval
t_start = 0;
t_stop = 4;
%Initial Condition
IC=.2;
%Solving the ODE
[t,y]=ode45(@(t,Y) rhs(t,Y,params), [t_start t_stop], IC); 
%plotting the ODE
plot(t,y,LineWidth=2)
xlabel('Time')
ylabel('Y')
%Define the right-hand-side function
function f=rhs(t,Y,params)
r=params(1);
k=params(2);
y=Y;
f=r*y*(k-y);
end