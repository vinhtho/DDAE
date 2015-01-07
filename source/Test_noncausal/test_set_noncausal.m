clear all; close all; clc;

%Exa_3_2_KunM06
Exa_noncausal_7

tspan = [0,10];

%option.MinStep = 0.01;
%option.MaxStep = 0.001;

options.StrIdx=1;
options.Shift=1;
%options.MaxStrIdx=1;

options.MaxShift=1;
options.Iter=1000;
options.MaxStrIdx = 3;

options.IsConst = 0;

sourcefolder = genpath('../');
addpath(sourcefolder);

tic
[t,x,info] = colddae(E,A,B,f,tau,phi,tspan,options);
toc

n = length(x(:,1));
X_e = xe(t);

figure(); clf;

for i=1:n
    subplot(2,n,i);
    plot(t,x(i,:))
    legend(strcat('x_',num2str(i)));

    subplot(2,n,n+i)
    plot(t,x(i,:)-X_e(i,:))
    legend(strcat('Abs. error of x_',num2str(i)));
end

rmpath(sourcefolder)    