addpath('C:\Users\Vinh Tho Ma\Documents\GitHub\DDAE')

tau = [1,1.5];

% exact solution
xe = @(t) [t; sin(t); exp(-t)];

% derivative of xe
xed = @(t) [1; cos(t); -exp(-t)];

% delayed xe
xet = @(t) [xe(t-tau(1)), xe(t-tau(2))];

F1 = @(t,x,xd,xt) xd(1) - x(1) - xt(1,:)*[1;1];
F2 = @(t,x,xt) x(2:3) + xt(2:3,:)*[1;2];

F1 = @(t,x,xd,xt) F1(t,x,xd,xt)-F1(t,xe(t),xed(t),xet(t));
F2 = @(t,x,xt) F2(t,x,xt)-F2(t,xe(t),xet(t));

F3 = 0;

% other parameters
t0 = 0;
x0 = xe(t0);
phi = @(t) xe(t);
h= 0.1;
N = 1000;
tol = 1e-7;

[t,x] = radar5(F1,F2,F3,t0,x0,tau,phi,h,N,tol);

semilogy(t,abs((x-xe(t))./xe(t)));