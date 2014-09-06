% DDAE TEST PROBLEM WITH ADDITIONAL DIFFERENTIAL EQUATIONS
% strangeness-index: 1
% shift index = 1

P=rand(2);
Q=rand(2);

E=@(t) P*[
    1   0   
    0   0   
    ]*Q;
A=@(t) P*[
    0   0  
    1   0  
    ]*Q;
B=@(t) P*[
    0   1   
    0   1   
    ]*Q;

phi=@(t)[
    t
    sin(t)
    ];

tau=1;
t0=0;
% consistent initial vector
x0=[0 0]';

% the exact solution
xe = @(t)[t;sin(t)];
xed = @(t)[1; cos(t)];

E0=E(0);
A0=A(0);
B0=B(0);
f=@(t) E0*[ 1; cos(t)] - A0*[t;sin(t)] - B0*[t-tau;sin(t-tau)] ;


options.Shift = 1;
options.StrIdx = 1;

options.MaxShift = 1;
options.MaxStrIdx = 1;

options.isConst = 1;


tic
% [t,x]=solve_shifted_lddae(E,A,B,f,tau,phi,[0,10],options);
[t,x]=solve_ddae({E,A,B,f,tau,phi},[0,10],options);
toc

semilogy(t,abs(x-xe(t)));
title('absolute error')