E=@(t)[1 0; 0 0];
A=@(t)[0 0; 1 0];
B=@(t)[0 1; 1 1];

phi = @(t) [sin(t);cos(t)];

tau = @(t) exp(-t);

f=@(t)E(t)*[cos(t);-sin(t)]-A(t)*phi(t)-B(t)*phi(t-tau(t));

tspan=[0,10];

options.StepSize=0.1;
options.StrIdx=1;
options.Shift=1;
options.MaxStrIdx=1;
options.MaxShift=1;

tic
[t,x] = solve_varshifted_lddae(E,A,B,f,tau,phi,tspan,options);
toc

close all
subplot(2,1,1)
semilogy(t,abs(x-phi(t)));
title('absolute error')
subplot(2,1,2)
plot(t,tau(t),tspan,[options.StepSize,options.StepSize])
legend('tau(t)','h')