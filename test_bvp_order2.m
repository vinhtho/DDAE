E = [
    0 0 1
    0 0 0
    0 1 0
    0 0 0
    ];

xe = @(t) [sin(t);cos(t);t/2];
xed = @(t) [cos(t);-sin(t);1/2];
phi = xe;
tau = 1;

%  b = @(t)(abs(t)<1).*exp(-1./(t.^2-1).^2);
%  a = @(t) b(t-3*tau);
% a = @(t) (t>2*tau).*(1-cos(t-2*tau));
a = @(t) t>tau;

A = @(t) [
    0    1 0 
    0    0 1
    0    0 0
    a(t) 0 0
    ];

B = [
    0 0 0
    1 0 0
    0 0 0
    0 0 0
    ];

f = @(t) E*xed(t)-A(t)*xe(t)-B*xe(t-tau);

E=@(t)E;
B=@(t)B;

t0 = 0;
t_end = 2 * tau;
h = 0.1;
N = floor((t_end-t0)/h);
tspan=[t0,t_end];


options.Shift = 1;
options.StrIdx = 1;
options.MaxStrIdx = 1;
options.MaxShift = 1;

tic
% [t,x] = solve_ddae({E,A,B,f,tau,phi},tspan,options);
[t,x] = solve_advanced_lddae(E,A,B,f,tau,phi,tspan,options);
% [t,x] = solve_varshifted_lddae(E,A,B,f,@(t)tau,phi,tspan,options);
toc

clf; figure(1);
subplot(2,2,1);
plot(t,x(1,:),t,x(2,:),'r-.','LineWidth',2);
legend('x_1','x_2')

err = abs(x-xe(t))./abs(xe(t));
subplot(2,2,2);
semilogy(t,err(1,:),t,err(2,:),'r-.','LineWidth',2);
legend('error x_1','error x_2')
