function test_nddae2(version)
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

%non-othogonal matrices lead to a larger errors
P=orth(rand(2));
Q=orth(rand(2));

E= P*[
    1   0   
    0   0   
    ]*Q;
A= P*[
    0   0  
    1   0  
    ]*Q;
B= P*[
    0   1   
    0   1  
    ]*Q;

tau=1;

% the exact solution
%
% cosine and sine
xe = @(t)[cos(t);sin(t)];
xed = @(t)[-sin(t); cos(t)];
xedd = @(t) -xe(t);
xeddd = @(t) -xed(t);
xedddd = xe;
%
% constant and linear function
% xe = @(t)[1*ones(size(t));t];
% xed = @(t)[zeros(size(t));ones(size(t))];
% xedd = @(t)zeros(2,numel(t));
% xeddd = xed;
% xedddd = xed;

f=@(t)  E*xed(t) - A*xe(t) - B*xe(t-tau);
df=@(t) E*xedd(t) - A*xed(t) - B*xed(t-tau);
ddf=@(t) E*xeddd(t) - A*xedd(t) - B*xedd(t-tau);
dddf=@(t) E*xedddd(t) - A*xeddd(t) - B*xeddd(t-tau);

F = { @(t,X,Xt) -E*X(3:4)+A*X(1:2)+B*Xt(1:2)+f(t)
      @(t,X,Xt) -E*X(5:6)+A*X(3:4)+B*Xt(3:4)+df(t)
      @(t,X,Xt) -E*X(7:8)+A*X(5:6)+B*Xt(5:6)+ddf(t)
      @(t,X,Xt) -E*X(9:10)+A*X(7:8)+B*Xt(7:8)+dddf(t)
      };


Phi = {xe,xed,xedd,xeddd,xedddd};

tspan = [0,5];

options.StrIdx = 1;
options.Shift = 1;
options.StepSize = 0.1;

% the dimension of x
n = 2;

tic
[t,x] = solve_nddae2(version,F,Phi,tau,tspan,n,options);
toc

close
subplot(2,1,1)
plot(t,x(1:2,:))
subplot(2,1,2)
semilogy(t,abs(x(1:2,:)-xe(t)))

disp(['maximal absolute error = ',num2str(max(max(abs(x(1:2,:)-xe(t)))))])