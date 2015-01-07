% System is noncausal of shift index 1
% strangeness index 2
% not of advanced type

P = @(t)[1 t t^2 t^3
         0 1 t   t^2
         0 0 t^2+1 t
         0 0 0     1];

E = [0 1 0
     0 0 1
     0 0 0
     0 1 1];

A=[1   0   0
   0   1   0
   0   0   0
   1   1   0];
    
B=[1   0   0
   0   0   0
   0   0   1
   1   0   1];

xe = @(t) [cos(t);sin(t);atan(t)];
xed= @(t) [-sin(t);cos(t);1/(1+t.^2)];

phi = xe;
tau = @(t) 1-0.5*sin(t);

f = @(t) P(t) * ( E*xed(t)-A*xe(t)-B*xe(t-tau(t)) );
E = @(t) P(t) * E;
A = @(t) P(t) * A;
B = @(t) P(t) * B;

tspan=[0,10];
