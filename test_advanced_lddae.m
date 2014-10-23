% DDAE TEST PROBLEM 14
% strangeness-index: 3

% advanced problem!!!

phi=@(t) [sin(t);cos(t);-sin(t);-cos(t)];
tau=1;
t0=0;
x0=phi(t0);

% exact solution
xe = @(t) phi(t);
xed = @(t) [cos(t);-sin(t);-cos(t);sin(t)];
xet = @(t) phi(t-tau);

E0=[
    0   1   0   0
    0   0   1   0
    0   0   0   1
    0   0   0   0
    0   1   1   1
    ];
A0=[
    1   0   0   0
    0   1   0   0
    0   0   1   0
    0   0   0   1
    1   1   1   1
    ];
B0=[
    0   0   0   0
    0   0   0   0
    0   1   0   0
    0   0   1   0
    0   1   1   0
    ];

E=@(t)E0;
A=@(t)A0;
B=@(t)B0;

% f=@(t) E(t)*xed(t)-A(t)*xe(t)-B(t)*xet(t);
f=@(t) E0*[cos(t);-sin(t);-cos(t);sin(t)] - A0*[sin(t);cos(t);-sin(t);-cos(t)] - B0*[sin(t-tau);cos(t-tau);-sin(t-tau);-cos(t-tau)];


options.StrIdx = 3;
options.MaxStrIdx = 10;
options.isConst = 1;

tic
[t,x]=solve_advanced_lddae(E,A,B,f,tau,phi,[0,3],options);
toc
semilogy(t,abs(x-xe(t))./abs(xe(t)));