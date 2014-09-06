
E = [
    1   0   0
    0   0   1
    0   0   0
    1   0   1
    ];

A = [
    1   0   0
    0   1   0
    0   0   1
    1   1   1
    ];

B = [
    1   0   0   2   0   0
    1   1   1   2   2   2
    0   0   1   0   0   0
    2   1   2   4   2   2
    ];


xe  = @(t) [sin(t); cos(t);atan(t)];
xed = @(t) [cos(t);-sin(t);1/(1+t.^2)];

phi = xe;
tau = [ 1 2 ];

f = @(t) E*xed(t)-A*xe(t)-B*[xe(t-tau(1));xe(t-tau(2))];
E = @(t) E;
A = @(t) A;
B = @(t) B;

tspan=[0,10];

options.StrIdx = 1;
options.StepSize = 0.2;

tic
% [t,x] = solve_lddae(E,A,B,f,tau,phi,tspan,options);
[t,x] = solve_ddae({E,A,B,f,tau,phi},tspan,options);
toc

semilogy(t,abs(x-xe(t)))
title('absolute error')