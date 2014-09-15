xe = @(t) [sin(t);cos(t)];
xed = @(t) [cos(t);-sin(t)];

E = [
    1 0
    0 0
    0 0
    ];

tau = 1;
phi = xe;

%  b = @(t)(abs(t)<1).*exp(-1./(t.^2-1).^2);
%  a = @(t) b(t-3*tau);
% a = @(t) (t>2*tau).*(1-cos(t-2*tau));
a = @(t) t>2*tau;

A = @(t) [
    0 0
    0 0
    0 a(t)
    ];

B = [
    0 0
    0 1
    0 0
    ];

f = @(t) E*xed(t)-A(t)*xe(t)-B*xe(t-tau);
E=@(t)E;
B=@(t)B;

tspan=[0,3*tau];

options.MaxStrIdx = 4;

tic
[t,x] = solve_ddae({E,A,B,f,tau,phi},tspan,options);
toc

semilogy(t,abs(x-xe(t)));