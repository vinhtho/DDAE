clear
E=@(t)[
    1 -1
    0 0];
A=@(t)[
    1 0
    0 1
    ];
B=@(t)[
    0 0
   -1 0
    ];
f=@(t)zeros(2,1);
tau=@(t)1;
phi=@(t)[1;1];
tspan=[0,4];

xe1 = @(t) exp(t);
xe2 = @(t) (t-1).*exp(t-1)+xe1(t);
xe3 = @(t) 0.5*(t.^2-2*t).*exp(t-2)+xe2(t);
xe4 = @(t) 1/6*(t.^3-3*t.^2-3*t+9).*exp(t-3)+xe3(t);
xe=@(t) xe1(t).*(0<=t).*(t<=1) + xe2(t).*(1<t).*(t<=2) + xe3(t).*(2<t).*(t<=3) + xe4(t).*(3<t).*(t<=4);

options.Iter=100;
options.MaxStrIdx=0;
[t,x,info] = colddae_causal(E,A,B,f,tau,phi,tspan,options);
semilogy(t,abs(x(1,:)-xe(t))./abs(xe(t)),'o-')
grid
xlabel('t')
legend('rel. error of x_1(t)','rel. error of x_2(t)')