function test_DDAE_with_add_diff_eq_01( )
%TEST_DDAE_WITH_ADD_DIFF_EQ_01 Summary of this function goes here
%   Detailed explanation goes here


P=eye(2);
Q=eye(2);

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

tspan = [0,2];

options.StrIdx = 3;
options.Shift = 1;
options.StepSize = 0.1;

% the dimension of x
n = 2;

[t,x] = solve_nddae2(F,Phi,tau,tspan,n,options);

err = max(max((abs(x(1:2,:)-xe(t)))));

assert(err<1e-2)

end

