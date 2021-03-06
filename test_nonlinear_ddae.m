% test radar5

clear
close all

F1=@(t,x,xd,xt) -xd(1) -sin(x(1)*x(2)) - [0 1 0 1]*xt;
F2=@(t,x,   xt)       x(2) - [1 0 1 0]*xt + sin(t); 

tau = @(t) [1-0.5*sin(t),t+1];
phi = @(t) [1;1];

options.NGrid=100;
options.tolR = 1e-7;
options.x0 = [1;2];

[t,x]=solve_nddae(F1,F2,tau,phi,[0,10],options);

plot(t,x)
title('solution')
