function [t,x]=solve_nddae2(version,F,Phi,tau,tspan,options)
%                                    .
% We want to solve the DDAE F(t,x(t),x(t),x(t-tau)) = 0 with constant delay
% tau>0.
% 
% The input parameter F is a struct containing the function F in the first 
% entry and its i-th time derivative at the (i+1)-st entry as a function
% handle (evaluation by feval(F{i},t,X,Xt)). Every function F{i} has t, X 
% and Xt as an input parameter, where X(j*n+1:(j+1)*n) contains the
% (j-1)-st derivative of x, same for Xt and x(t-tau).
%
% Analogously Phi is a struct containing the i-th derivative of the history
% function at the (i+1)-st entry having only t as an input.
%
%   Use version to switch between several systems.
%   
%   version=1: Inflated system without any additional equations
% 
%   version=2: Inflated system and F(t,x,Dh x) = 0
%                                  .
%   version=3: Inflated system and x = Dh x 
%                                                     .
%   version=4: Inflated sytem and F(t,x,Dh x) = 0 and x = Dh x 
% 


% the shift and strangeness index have to be known
K = 0;
KMax = 3;
mu = 0;
muMax = 3;
N = 99;
h = diff(tspan)/N;
tolA = 1e-7; % not used right now
tolR = 1e-7;
if exist('options','var')
    if isfield(options,'isConst')   isConst=options.isConst; end
    if isfield(options,'AbsTol')    tolA=options.AbsTol; end
    if isfield(options,'MaxShift')  KMax=options.MaxShift; end
    if isfield(options,'MaxStrIdx') muMax=options.MaxStrIdx; end
    if isfield(options,'NGrid')     N=options.NGrid-1; h=diff(tspan)/N; end
    if isfield(options,'RelTol')    tolR=options.RelTol; end
    if isfield(options,'Shift')     K=max(0,options.Shift); KMax = max(KMax,K); end
    if isfield(options,'StepSize')  h=options.StepSize; N=floor(diff(tspan)/h); end
    if isfield(options,'StrIdx')    mu=options.StrIdx; muMax = max(muMax,mu); end
    if isfield(options,'x0')        x0=options.x0; end
else
    options = {};
end

% the dimension
n=length(feval(Phi{1},0));

% the output parameters
x = zeros((mu+2)*n*(K+1),N+1);
t = tspan(1):h:tspan(2);

% Actually, we are cheating here, because we don't know x(t) for t>t0.
% However, in most of our test examples x(t) = phi(t) for all t.
for i = 0:K
    for j = 0:mu
        x(j*n+1+i*(mu+2)*n:(j+1)*n+i*(mu+2)*n,1) = feval(Phi{j+1},tspan(1)+i*tau);
    end
end


% main loop
for k = 1:N
    %                          .            (mu)
    % Computing Xt = [x(t-tau);x(t-tau);...;x (t-tau)]
    Xt = zeros((mu+1)*n,1);
    if t(k)+h-tau<tspan(1)
        for j = 0:mu
            Xt(j*n+1:(j+1)*n,1) = feval(Phi{j+1},t(k)+h-tau);
        end
    else
        k_tau = find(t(k)+h-tau<=t,1);
        % we use linear interpolation between the x-values corresponding to
        % t(k_tau) and t(k_tau+1)
        Xt(1:(mu+1)*n,1) = x(1:(mu+1)*n,k_tau)+(x(1:(mu+1)*n,k_tau+1)-x(1:(mu+1)*n,k_tau))*(t(k)+h-tau-t(k_tau))/(t(k_tau+1)-t(k_tau));
    end
    
    % 
    bdf_order = 6;
    % The matrix xx contains the bdf_order x-values to x(:,k+1),
    % i.e. x(:,k+1-bdf_order), x(:,k+1-bdf_order-1), x(:,k+1-bdf_order-2),
    % ..., x(:,k). For time points smaller than t0 we use the history
    % function.
    if k >= bdf_order+1
        xx = x(:,k-bdf_order+1:k);
    else
        xx = zeros((mu+2)*n*(K+1),bdf_order);
        for l = 1:bdf_order
            if t(k)-(l-1)*h<t(1)
                for j = 0:K
                    xx(j*n+1:(j+1)*n,bdf_order-l+1) = feval(Phi{1},t(k)-(l-1)*h);
                end
            else
                xx(:,bdf_order-l+1) = x(:,k-l+1);
            end
        end
    end
    
    % xdata is just a dummy variable for lsqcurvefit
    FF = @(X,xdata)[
        form_big_system(F,t(k)+h,X,Xt,K,mu,n,tau);
        form_bdf_system(version,F,t(k)+h,xx,X,Xt,h,n,bdf_order)
        ];
    
    switch version
        case 1
            y = zeros((K+1)*(mu+1)*n,1);
        case 2
            y = zeros((K+1)*(mu+1)*n+n,1);
        case 3
            y = zeros((K+1)*(mu+1)*n+n,1);
        case 4
            y = zeros((K+1)*(mu+1)*n+2*n,1);
    end
    
    % Now solve FF(X)=0
    warning off
    x(:,k+1) = lsqcurvefit(FF,x(:,k),{},y,[],[],optimset('Display','Off'));
    warning on
end

function FF1 = form_big_system(F,t,X,Xt,K,mu,n,tau)
% X  has length (mu+2)*n*(K+1)
% Xt has length (mu+1)*n
FF1 = zeros((mu+1)*n*(K+1),1);
for j = 0:mu
    FF1(j*n+1:(j+1)*n) = feval(F{j+1},t,X(1:(j+2)*n),Xt(1:(j+1)*n));
end
for i = 1:K
    for j = 0:mu
        FF1(j*n+1+i*(mu+1)*n:(j+1)*n+i*(mu+1)*n) = feval(F{j+1},t+i*tau,X((1:(j+2)*n)+i*(mu+2)*n),X((1:(j+1)*n)+(i-1)*(mu+2)*n));
    end
end

function FF2 = form_bdf_system(version,F,t,xx,X,Xt,h,n,bdf_order)
switch bdf_order
    case 1
        Dh_X = (X-xx)/h;
    case 2
        Dh_X = (3/2*X+xx*[1/2; -2])/h;
    case 3
        Dh_X = (11/6*X+xx*[-1/3; 3/2; -3])/h;
    case 4
        Dh_X = (25/12*X+xx*[1/4; -4/3; 3; -4])/h;
    case 5
        Dh_X = (137/60*X+xx*[-1/5; 5/4; -10/3; 5; -5])/h;
    case 6
        Dh_X = (49/20*X+xx*[1/6; -6/5; 15/4; -20/3; 15/2; -6])/h;
end

switch version
    case 1
        FF2 = [];
    case 2
        FF2 = feval(F{1},t,[X(1:n);Dh_X(1:n)],Xt);
    case 3
        FF2 = X(n+1:2*n)-Dh_X(1:n);
    case 4
        FF2 = zeros(2*n,1);
        FF2(1:n) = feval(F{1},t,[X(1:n);Dh_X(1:n)],Xt);
        FF2(n+1:2*n) = X(n+1:2*n)-Dh_X(1:n);
end
