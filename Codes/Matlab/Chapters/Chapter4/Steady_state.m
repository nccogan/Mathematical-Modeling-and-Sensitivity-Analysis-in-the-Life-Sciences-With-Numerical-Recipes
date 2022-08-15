clear
close all
params=[.2,.1,.2,.1];
%plot contours
 [XX,YY]=meshgrid([-2:.1:2],[-2:.1:2]);
 [m,n]=size(XX);
 for i=1:m
     for j=1:n
         temp=rhs([XX(i,j),YY(i,j)],params);
         A(i,j)=temp(1);
         B(i,j)=temp(2);
     end
 end
 contour(XX,YY,A,[-.0001,.0001],'k','Linewidth',2)
 hold on
 contour(XX,YY,B,[-.0001,.0001],'c','Linewidth',2)
%=======================================================
%find Steady-States: Using contour to note there are two
%=======================================================
solution1 = fsolve(@(Y)rhs(Y,params), [10,.1])
solution2= fsolve(@(Y)rhs(Y,params), [-10,-.1])
plot(solution1(1),solution1(2),'x','MarkerSize',15,'LineWidth',2)
plot(solution2(1),solution2(2),'x','MarkerSize',15,'LineWidth',2)

%====================================
%Define the right hand side functions
%=====================================
function f=rhs_x(Y,params)
alpha=params(1);
beta=params(2);
delta=params(3);
gamma=params(4);
x=Y(1);
y=Y(2);
f=alpha.*x.*y-delta*x;
end

function g=rhs_y(Y,params)
alpha=params(1);
beta=params(2);
delta=params(3);
gamma=params(4);
x=Y(1);
y=Y(2);
g =-beta.*x.*y+gamma*y;
end

function k=rhs(Y,params)
alpha=params(1);
beta=params(2);
delta=params(3);
gamma=params(4);
x=Y(1);
y=Y(2);
k = [rhs_x(Y,params),rhs_y(Y,params)];
end

