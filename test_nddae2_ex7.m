% test_nddae2_ex7
%
%   Use version to switch between several systems.
%   
%   version=1: Inflated system without any additional equations.
%
%   version=2: Inflated system and F(t,x,Dh x) = 0.
%                                  .
%   version=3: Inflated system and x = Dh x 
%                                                     .
%   version=4: Inflated sytem and F(t,x,Dh x) = 0 and x = Dh x 
%

tau=1;

% the exact solution
%
xe = @(t)[1-exp(t);t;t-1];
xed = @(t)[-exp(t);1;1];
xedd = @(t) [-exp(t);0;0];

F = { @(t,x,xt) [
        x(6)-xt(1)-exp(t-1)
%         x(5)-x(1)-exp(t)
        x(4)*(x(5)-x(1)-exp(t))+x(2)-t
        x(3)-xt(2)
        ],
      @(t,x,xt) [
        x(9)-xt(4)-exp(t-1)
        x(7)*(x(5)-x(1)-exp(t))+x(4)*(x(8)-x(4)-exp(t))+x(5)-1
        x(6)-xt(5)
        ]
     };
 
F1 = @(t,x,xd,xt) [
        xd(3)-xt(1)-exp(t)
        xd(1)*(xd(2)-x(1)-exp(t))+x(2)-t
        ];
F2 = @(t,x,xt) x(3)-xt(2);


Phi = {xe,xed};

tspan = [0,5];

options.StrIdx = 1;
options.Shift = 1;
options.StepSize = 0.1;


version = 1;
n = 3;

tic
[t,x] = solve_nddae2(version,F,Phi,tau,tspan,n,options);
toc

close
subplot(2,1,1)
plot(t,x(1:3,:))
subplot(2,1,2)
semilogy(t,abs(x(1:3,:)-xe(t))./abs(xe(t)))