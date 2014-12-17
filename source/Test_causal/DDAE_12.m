% DDAE TEST PROBLEM 12
% strangeness-index: 2
% 2 constant delay

E=@(t) [
    0   1   0
    0   0   1
    0   0   0
    ];

A=@(t) [
    1   0   0
    0   1   0
    0   0   1
    ];

B1=@(t) [
    0   1   0
    0   0   0
    0   0   0
    ];

B2=@(t) [
    0   0   1
    0   0   0
    0   0   0
    ];

B=@(t)[B1(t),B2(t)];

phi=@(t)[
    exp(t)
    t
    sin(t)
    ];

tau=@(t)[pi/100, pi/2];
t0=0;
t_end=10;

h=0.01;
N =1000;
% inconsistent initial vector
%x0=[-3 0 -1]';
x0=[1 0 0]';

f=@(t) E(t)* [exp(t) 1 cos(t)]' - A(t) * phi(t) - B1(t)*phi(t-pi/100) - B2(t)*phi(t-pi/2);
eps0=0.01;

% the exact solution
xe = @(t)phi(t);