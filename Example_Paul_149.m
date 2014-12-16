clear;
m=2;
E=@(t)eye(4);
A=@(t)[
    0 0 1 0
    0 0 0 1
    0 -2*m 0 0
    -2*m 0 0 0
    ];
B=@(t)[
    0 0 0 0
    0 0 0 0
    (1+m^2)*(-1)^m 0 0 0
    0 (1+m^2)*(-1)^m 0 0
    ];
f=@(t)zeros(4,1);
tau=@(t)pi;
phi=@(t)[
    sin(t).*cos(m*t)
    cos(t).*sin(m*t)
    cos(t).*cos(m*t)-m*sin(t).*sin(m*t)
    m*cos(t).*cos(m*t)-sin(t).*sin(m*t)
    ];
phi12=@(t)[
    sin(t).*cos(m*t)
    cos(t).*sin(m*t)
    ];
tspan=[0,5];
options.Iter=100;
options.MaxStrIdx=0;
[t,x,info] = colddae_causal(E,A,B,f,tau,phi,tspan,options);
semilogy(t,abs(x(1:2,:)-phi12(t))./abs(phi12(t)),'o-')
grid
xlabel('t')
legend('rel. error of x_1(t)','rel. error of x_2(t)')