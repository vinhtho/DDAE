E=@(t)[1 0; 0 0; 1 0];
A=@(t)[0 0; 1 0; 1 0];
B=@(t)[0 1; 1 1; 1 2];

% phi = @(t) [atan(t);-atan(t)];
% dphi = @(t) [1/(1+t.^2);-1/(1+t.^2)];

phi = @(t)[sin(t);cos(t)];
dphi= @(t)[cos(t);-sin(t)];

tau = @(t) exp(-t);
% tau=@(t)sin(t)+2;

f=@(t)E(t)*dphi(t)-A(t)*phi(t)-B(t)*phi(t-tau(t));

tspan=[0,10];

options.InitStep=0.1;
% options.MaxStep=0.5;
options.MinStep=0.1;

options.StrIdx=1;
options.Shift=1;
options.MaxStrIdx=1;
options.MaxShift=1;

options.RelTol=1e-3;
options.AbsTol=1e-3;
options.LagTol=1e-2;

tic
[t,x,info] = solveCausalLDDAE(E,A,B,f,tau,phi,tspan,options);
toc

close all
subplot(2,1,1)
semilogy(t,abs(x-phi(t)),'o-');
title('absolute error')
subplot(2,1,2)
plot(t,tau(t),t(1:end-1),diff(t),'o-')
legend('tau(t)','h')