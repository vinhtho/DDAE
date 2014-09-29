% F = { @(t,X,Xt)[
%                     X(4)-X(1)
%                     X(3)*(X(4)-X(1))+Xt(2)
%     ];@(t,X,Xt)[
%                     X(6)-X(3)
%                     X(5)*(X(4)-X(1))+X(3)*(X(6)-X(3))+Xt(4)
%     ]};
% 
% Phi = {@(t)[sin(t),cos(t)],@(t)[cos(t),-sin(t)]};

P=rand(2);
Q=rand(2);

E= P*[
    1   0   
    0   0   
    ]*Q;
A= P*[
    0   1  
    1   0  
    ]*Q;
B= 0*P*[
    0   1   
    0   0   
    ]*Q;

tau=1;

f=@(t)  E*[ 1; cos(t)] - A*[t;sin(t)] - B*[t-tau;sin(t-tau)];
df=@(t) E*[ 0;-sin(t)] - A*[1;cos(t)] - B*[1;cos(t-tau)];


F = { @(t,X,Xt) -E*X(3:4)+A*X(1:2)+B*Xt(1:2)+f(t),
      @(t,X,Xt) -E*X(5:6)+A*X(3:4)+B*Xt(3:4)+df(t)
      };


Phi = {@(t)[t,sin(t)],@(t)[1,cos(t)]};


% the exact solution
xe = @(t)[t;sin(t)];
xed = @(t)[1; cos(t)];

tspan = [0,10];

options.StrIdx = 0;
options.Shift = 0;

[t,x] = solve_nddae2(F,Phi,tau,tspan,2,options);

plot(t,x(1:2,:))