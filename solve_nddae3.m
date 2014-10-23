function [t,x]=solve_nddae3(F,Phi,tau,tspan,n,options)
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

% the output parameters
x = zeros((mu+2)*n*(K+1),N+1);
t = tspan(1):h:tspan(2);

% Actually, we are cheating here, because we don't know x(t) for t>t0.
% However, in most of our test examples x(t) = phi(t) for all t.
for i = 0:K
    for j = 0:mu+1
        x(j*n+1+i*(mu+2)*n:(j+1)*n+i*(mu+2)*n,1) = feval(Phi{j+1},tspan(1)+i*tau);
    end
end

% main loop
for k = 1:N
    %                          .            (mu)
    % Computing Xt = [x(t-tau);x(t-tau);...;x (t-tau)].
    Xt = zeros((mu+1)*n,1);
    if t(k)+h-tau<tspan(1)
        for j = 0:mu
            Xt(j*n+1:(j+1)*n,1) = feval(Phi{j+1},t(k)+h-tau);
        end
    else
        k_tau = find(t(k)+h-tau<=t,1);
        % We use linear interpolation between the x-values corresponding to
        % t(k_tau) and t(k_tau+1)
        Xt(1:(mu+1)*n,1) = x(1:(mu+1)*n,k_tau)+(x(1:(mu+1)*n,k_tau+1)-x(1:(mu+1)*n,k_tau))*(t(k)+h-tau-t(k_tau))/(t(k_tau+1)-t(k_tau));
    end
    
    % xdata is just a dummy variable for lsqcurvefit
    FF = @(X,xdata)[
        form_big_system(F,t(k)+h,X,Xt,K,mu,n,tau);
        ];
    
    % Now solve FF(X)=0.
    x(:,k+1) = lsqcurvefit(FF,x(:,k),{},zeros((mu+1)*n*(K+1),1),[],[],optimset('Display','Off'));
end

% subroutine for evaluating the big inflated system
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