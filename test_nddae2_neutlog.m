%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    3D Gause-type predator-pray model
%
% On neutral delay logistic {G}ause-type predator-prey systems
% Y. Kuang
% Dynam. Stability Systems, vol. 6, pp. 173--189, 1991
%
%
%   y1�(t) = y1(t) (1 - y1(t-tau ) - rho y3(t-tau ) - y2(t) F(y1(t))
%
%   y2�(t) = y2(t) (F(y1(t)) - alpha)
%
%        0 = y1(t) (1 - y1(t-tau) - rho y3(t-tau)) - y2(t) F(y1(t)) - y3(t)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all
clear all
clc

% SETTING THE PARAMETERS
alpha   = 0.1;
rho     = 2.9;
tau     = 0.42;

% THE DDAE
f       = @(x) x^2/(x^2+1);
df      = @(x) (2*x*(x^2+1)-2*x^3)/(x^2+1)^2;

F = {
    @(t,x,xt) [
    -x(4)+x(1)*(1-xt(1)-rho*xt(3))-x(2)*f(x(1))
    -x(5)+x(2)*(f(x(1))-alpha)
    x(1)*(1-xt(1)-rho*xt(3))-x(2)*f(x(1))-x(3)
    ],
    @(t,x,xt) [
    -x(7)+x(4)*(1-xt(1)-rho*xt(3))+x(1)*(-xt(4)-rho*xt(6))-x(5)*f(x(1))-x(2)*df(x(1))*x(4)
    -x(8)+x(5)*(f(x(1))-alpha)+x(2)*(df(x(1))*x(4))
    x(4)*(1-xt(1)-rho*xt(3))+x(1)*(-xt(4)-rho*xt(6))-x(5)*f(x(1))-x(2)*df(x(1))*x(4)-x(6)
    ]};

Phi     = {
    @(t) [
    0.33-t/10
    2.22+t/10
    -1/10
    ]
    @(t)zeros(3,1);
    };

% the appropriate function handles for solve_nddae (radar5)
F1 = @(t,x,xd,xt) [
                    -xd(1)+x(1)*(1-xt(1)-rho*xt(3))-x(2)*f(x(1))
                    -xd(2)+x(2)*(f(x(1))-alpha) ];
F2 = @(t,x,xt)             x(1)*(1-xt(1)-rho*xt(3))-x(2)*f(x(1))-x(3);
phi = Phi{1};


tspan = [0,6];

 x0t= feval(Phi{1},tspan(1)-tau);
 x01= 0.33;
 x02= 0.22;
% %x03= -0.1;
% %x02= -(x03-x01*(1-x0t(1)-rho*x0t(3)))/(f(x01));
 x03= x01*(1-x0t(1)-rho*x0t(3))-x02*f(x01);
% options.x0 = [x01;x02;x03];

options.Shift = 0;
options.StrIdx = 0;
options.StepSize = 0.1;

version = 3;

n=3;

disp('COMMENCING SOLVING PROCESS...')
tic
[t,x]=solve_nddae2(version,F,Phi,tau,tspan,options);
toc
disp('SOLVING COMPLETE.')
figure
plot(t,x(1:3,:))
legend('x_1(t)','x_2(t)','x_3(t)')