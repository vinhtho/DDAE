clear all; close all; clc;

DDAE_15;
% trouble: 1,2,4,5,6,7,9,11,12,13,14,14b,15

tspan = [0,10];
options.Iter=100;
options.MaxStrIdx = 3;

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
    