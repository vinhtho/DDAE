E = [0 1 0; 0 0 1; 0 0 0];
A = eye(3);
B = eye(3);

xe = @(t) [cos(t);sin(t);atan(t)];
xed= @(t) [-sin(t);cos(t);1/(1+t.^2)];

phi = xe;
tau = @(t) 1-0.5*sin(t);

f = @(t) E*xed(t)-A*xe(t)-B*xe(t-tau(t));
E = @(t) E;
A = @(t) A;
B = @(t) B;

tspan=[0,10];

options.StrIdx = 2;

tic
% [t,x] = solve_timeVarDelayed_strange_lddae(E,A,B,f,tau,phi,tspan,options);
[t,x] = solve_ddae({E,A,B,f,tau,phi},tspan,options);
toc

semilogy(t,abs(x-xe(t)))
title('absolute error')