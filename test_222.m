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
% the exact solution
xe = @(t)phi(t);

tau=@(t)[1, t/2+1];

t0=0;
t_end=10;

h=0.01;
N =1000;
% inconsistent initial vector
%x0=[-3 0 -1]';
x0=[1 0 0]';

f=@(t) E(t)* [exp(t) 1 cos(t)]' - A(t) * phi(t) - B1(t)*phi(t-1) - B2(t)*phi(t-t/2-1);
eps0=0.01;

F1=@(t,x,xdot,xtau) eye(2,3) * (-E(t)*xdot + A(t)*x + B1(t)*xtau(1:3) + B2(t)*xtau(4:6) + f(t));
F2=@(t,x,xtau) [0 0 1] * (A(t)*x + B1(t)*xtau(1:3) + B2(t)*xtau(4:6) + f(t));

options.NGrid=100;
options.tolR = 1e-7;
options.x0=x0;

[t,x]=solve_nddae(F1,F2,tau,phi,[0,10],options);

clf; figure(1);
subplot(1,2,1)
title('solution')
plot(t,x)

subplot(1,2,2)
title('absolute solution')
plot(t,abs(x-xe(t)))