% Over determined system
% shift index 1 if t < 2tau, 0 if t > 2tau
% strangeness index always 0
% of retarded type

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

%f = @(t) E*xed(t)-A(t)*xe(t)-B*xe(t-tau);
f = @(t)[cos(t);-cos(t-tau);-a(t)*cos(t)];

E=@(t)E;
B=@(t)B;

tspan=[0,10*tau];

options.MaxStrIdx = 0;

tic
% shifted solver, Radau IIA collocation works fine
% [t,x] = solve_ddae({E,A,B,f,tau,phi},tspan,options);

% BVP solver for noncausal does not work
[t,x] = solve_advanced_lddae(E,A,B,f,tau,phi,tspan,options);
toc

clf; figure(1);
subplot(2,2,1);
plot(t,x(1,:),t,x(2,:),'r-.','LineWidth',2);
legend('x_1','x_2')

err = abs(x-xe(t))./abs(xe(t));
subplot(2,2,2);
semilogy(t,err(1,:),t,err(2,:),'r-.','LineWidth',2);
legend('error x_1','error x_2')

% print -depsc bvp_noncausal.eps
% ! epstopdf bvp_noncausal.eps
% ! rm *.eps