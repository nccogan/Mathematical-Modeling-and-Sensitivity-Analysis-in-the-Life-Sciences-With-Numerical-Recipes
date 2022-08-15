%================================================================
% Define the QoI function
%================================================================
function SA = QoI(params)
    V0=6;
    m0=.053;
    h0=.595;
    n0=.317;
    tspan = linspace(0, 20, 500);
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
 %====================
     [t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)],[V0,m0,h0,n0]); 
     SA = max(y_solution(:,1));

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

end
