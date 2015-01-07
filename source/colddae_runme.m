% COLDDAE_RUNME: A script file for demonstrating the use of the DDAE solver
%                COLDDAE
clear,close all,clc

disp('Example 1: Multiple delays and strangeness-index 1 and default parameters.')
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
    0   0   0   0   0   0
    2   1   2   4   2   2
    ];
xe  = @(t) [sin(t); cos(t);atan(t)];
xed = @(t) [cos(t);-sin(t);1/(1+t.^2)];
phi = xe;
tau = @(t)[1 2];
f = @(t) E*xed(t)-A*xe(t)-B*[xe(t-1); xe(t-2)];
E = @(t) E;
A = @(t) A;
B = @(t) B;
tspan=[0,2*pi];
[t,x,info] = colddae(E,A,B,f,tau,phi,tspan);
info
semilogy(t,abs(x-xe(t)),'o-')
grid
title('Example 1: Absolute error.')
disp('The error becomes quite large at the end.')

disp('Press enter to continue.')
pause

disp('Example 2: Same system, but rel and abs tolerance are now 1e-7.')
options.RelTol=1e-7;
options.AbsTol=1e-7;
[t,x,info] = colddae(E,A,B,f,tau,phi,tspan,options);
info
semilogy(t,abs(x-xe(t))./abs(xe(t)),'o-')
grid
title('Example 1: Relative error')
disp('The error looks a little better now.')

disp('Press enter to continue.')
pause

