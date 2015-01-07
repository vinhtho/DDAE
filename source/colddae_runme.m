% COLDDAE_RUNME: A script file for demonstrating the use of the DDAE solver
%                COLDDAE
clear,close all,clc

% define the system
E0 = [
    1   0   0
    0   0   1
    0   0   0
    1   0   1
    ];
A0 = [
    1   0   0
    0   1   0
    0   0   1
    1   1   1
    ];
B0 = [
    1   0   0   2   0   0
    1   1   1   2   2   2
    0   0   0   0   0   0
    2   1   2   4   2   2
    ];
xe  = @(t) [sin(t); cos(t);atan(t)];
xed = @(t) [cos(t);-sin(t);1./(1+t.^2)];
phi = xe;
tau = @(t)[1 2];
f = @(t) E0*xed(t)-A0*xe(t)-B0*[xe(t-1); xe(t-2)];
E = @(t) E0;
A = @(t) A0;
B = @(t) B0;
tspan=[0,10];

disp('Example 1: Causal overdetermined system with multiple constant delays and') 
disp('strangeness-index 1 and default parameters.')
[t,x,info] = colddae(E,A,B,f,tau,phi,tspan);
info
subplot(2,1,1)
plot(t,x,'o',t,xe(t))
subplot(2,1,2)
semilogy(t,abs(x-xe(t)),'o-')
grid
title('solution and abs error')
disp('The error becomes quite large at the end.')

disp('Press enter to continue.')
disp(' ')
pause

disp('Example 2: Same as Example 1, but rel and abs tolerance are now 1e-9.') 
disp('(This might take some time.)')
options.RelTol=1e-9;
options.AbsTol=1e-9;
[t,x,info] = colddae(E,A,B,f,tau,phi,tspan,options);
info
subplot(2,1,1)
plot(t,x,'o',t,xe(t))
subplot(2,1,2)
semilogy(t,abs(x-xe(t)),'o-')
grid
disp('The error looks a little better now.')

disp('Press enter to continue.')
disp(' ')
pause

disp('Example 3: Same as Example 2, but since E,A,B,tau are constant,') 
disp('we can save some time by computing the strangenes-free form only once.')
options.IsConst = 1;
[t,x,info] = colddae(E,A,B,f,tau,phi,tspan,options);
info
subplot(2,1,1)
plot(t,x,'o',t,xe(t))
subplot(2,1,2)
semilogy(t,abs(x-xe(t)),'o-')
grid
disp('The computation should have taken less time now.')

disp('Press enter to continue.')
disp(' ')
pause

disp('Example 4: Same as Example 3, but now the derivative array is provided.')

M = @(t)[-A0,E0,zeros(4,3);zeros(4,3),-A0,E0];
P = @(t)[B0,zeros(4,6);zeros(4,6),B0];
xedd = @(t) [-sin(t);-cos(t);-2*t./(1+t.^2).^2];
g = @(t) [
    E0*xed(t)-A0*xe(t)-B0*[xe(t-1); xe(t-2)]
    E0*xedd(t)-A0*xed(t)-B0*[xed(t-1); xed(t-2)]
    ];

options.DArray = {M,P,g};

[t,x,info] = colddae(E,A,B,f,tau,phi,tspan,options);
info
subplot(2,1,1)
plot(t,x,'o',t,xe(t))
subplot(2,1,2)
semilogy(t,abs(x-xe(t)),'o-')
grid
disp('The error should be smaller, since we use exact derivatives for the regularization.')

disp('Press enter to continue.')
disp(' ')
pause

% new system
clear

E0 = [0 1 0
     0 0 1
     0 0 0
     0 1 1];
A0=[1   0   0
   0   1   0
   0   0   0
   1   1   0];
B0=[1   0   0
   0   0   0
   0   0   1
   1   0   1];
xe = @(t) [cos(t);sin(t);atan(t)];
xed= @(t) [-sin(t);cos(t);1/(1+t.^2)];
phi = xe;
tau = @(t) 1-0.5*sin(t);
f = @(t) E0*xed(t)-A0*xe(t)-B0*xe(t-tau(t));
E = @(t) E0;
A = @(t) A0;
B = @(t) B0;
tspan=[0,10];

disp('Example 5: Noncausal overdetermined system with single variable delay and') 
disp('strangeness-index 2, shift index 1 and default parameters.')
[t,x,info] = colddae(E,A,B,f,tau,phi,tspan);
info
subplot(2,1,1)
plot(t,x,'o',t,xe(t))
subplot(2,1,2)
semilogy(t,abs(x-xe(t)),'o-')
grid
title('solution and abs error')

disp('Press enter to continue.')
disp(' ')
pause

disp('Example 6: Same as Example 5, but strangeness- and shift-index are now set.')
clear options
options.StrIdx=2;
options.MaxStrIdx=2;
options.Shift=1;
options.MaxShift=1;
[t,x,info] = colddae(E,A,B,f,tau,phi,tspan,options);
info
subplot(2,1,1)
plot(t,x,'o',t,xe(t))
subplot(2,1,2)
semilogy(t,abs(x-xe(t)),'o-')
grid
title('solution and abs error')
disp('The computation should have taken less time now.')

disp('Press enter to continue.')
disp(' ')
pause


disp('Example 7: Same as Example 6, but Strangeness- and shift-index are now set')
disp('and derivative array is provided')

M = @(t)[-A0,E0,zeros(4,6);zeros(4,3),-A0,E0,zeros(4,3);zeros(4,6),-A0,E0];
P = @(t)[
    B0,zeros(4,6);
    zeros(4,3),B0*(1+0.5*cos(t)),zeros(4,3)
    zeros(4,3),-B0*0.5*sin(t),B0*(1+0.5*cos(t))^2];
xedd = @(t) [-sin(t);-cos(t);-2*t./(1+t.^2).^2];
xeddd= @(t) [-cos(t);sin(t);(6*t^2-2)/(t^2+1)^3];
g = @(t) [
    E0*xed(t)-A0*xe(t)-B0*xe(t-1+0.5*sin(t))
    E0*xedd(t)-A0*xed(t)-B0*xed(t-1+0.5*sin(t))*(1+0.5*cos(t))
    E0*xeddd(t)-A0*xedd(t)-B0*xedd(t-1+0.5*sin(t))*(1+0.5*cos(t))^2 + B0*xed(t-1+0.5*sin(t))*0.5*sin(t)
    ];

options.DArray = {M,P,g};

[t,x,info] = colddae(E,A,B,f,tau,phi,tspan,options);
info
subplot(2,1,1)
plot(t,x,'o',t,xe(t))
subplot(2,1,2)
semilogy(t,abs(x-xe(t)),'o-')
grid
title('solution and abs error')
disp('The computation should have taken less time now.')