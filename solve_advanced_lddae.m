function [t,x] = solve_advanced_lddae(E,A,B,f,tau,phi,tspan,options)
%solve_advanced_lddae
%
%   Solver for linear delay differential algebraic equations (DDAEs) with 
%   variable coefficients and single, constant delay tau > 0, i.e.
%           .
%       E(t)x(t) = A(t)x(t) + B(t)x(t-tau) + f(t),     t0 <= t <= tf,
%           x(t) = phi(t),                             t <= t0,
%   
%   based reformulating the DDAE as a boundary value problem and using 
%   Radau IIA collocation.
%                          .
%   The corresponding DAE Ex = Ax can have strangeness index bigger than 
%   zero (differentiation index bigger than one). However, note that bigger
%   errors might occur if the strangeness index is too big, because
%   derivatives are approximated by finite differences.
%
%   This solver can handle certain advanced DDAEs, i.e. the solution of
%   the DDAE x at the time t depends on derivatives of x(t-tau). Similar to
%   above, the higher the order of these derivatives is, the larger are the
%   errors. As the solution propagates in time, the order of the
%   derivatives can grow. In this case one must restrict the size of the
%   time interval.
%
%   INPUT   
%   -----   
%   E       fcn_handle      m-by-n leading matrix function
%   A       fcn_handle      m-by-n matrix function
%   B       fcn_handle      m-by-n matrix function for the delayed part
%   f       fcn_handle      m-by-1 inhomogeneity function
%   tau     double(1,1)     constant delay
%   phi     fcn_handle      history function
%   tspan   double(1,2)     the time interval
%   options struct          contains several fields with respective values,
%                           see options below
%
%   FIELDS FOR OPTIONS
%   ------------------
%   AbsTol      double      absolute tolerance
%   isConst     boolean     TRUE if E,A,B are constant, FALSE otherwise
%   MaxStrIdx   double      upper bound for the strangeness index
%   NGrid       double      no. of grid points for time discretization
%   RelTol      double      relative tolerance
%   StepSize    double      (constant) step size for the one-step method
%   StrIdx      double      strangeness index (if known)
%   x0          double(n,1) initial value (usually x0 = phi(t0))
%
%   OUTPUT      
%   ------  
%   t       double(1,Ngrid)   discretization of tspan
%   x       double(n,Ngrid)   approximate solution at the time nodes in t
%
%   References:
%   http://www.math.uni-bremen.de/zetem/cms/media.php/262/report9908.pdf
%   TODO more references

% set the options
isConst = 0;
mu0 = 0;
muMax = 3;
N = 99;
h = diff(tspan)/N;
tolA = 1e-7;
tolR = 1e-7;
x0 = phi(tspan(1));
if exist('options','var')
    if isfield(options,'isConst')   isConst=options.isConst; end
    if isfield(options,'AbsTol')    tolA=options.AbsTol; end
    if isfield(options,'MaxStrIdx') muMax=options.MaxStrIdx; end
    if isfield(options,'NGrid')     N=options.NGrid-1; h=diff(tspan)/N; end
    if isfield(options,'RelTol')    tolR=options.RelTol; end
    if isfield(options,'StepSize')  h=options.StepSize; N=floor(diff(tspan)/h); end
    if isfield(options,'StrIdx')    mu0=options.StrIdx; muMax = max(muMax,mu0); end
    if isfield(options,'x0')        x0=options.x0; end
else
    options = {};
end

% the initial time
t0 = tspan(1);

% the number of intervals of the form [(i-1)*tau,i*tau] contained in tspan
L = floor(diff(tspan)/tau);
N = floor(N/L);
options_new = options;
options_new.NGrid = N+1;

% the dimension
[m,n] = size(E(t0));

% the BVP is given by the following 
EE = @(s) E_reformulated(E,m,n,L,s,tau);
AA = @(s) A_reformulated(A,B,m,n,L,s,tau);
ff = @(s) f_reformulated(f,B,phi,m,n,L,s,tau);

% find a consistent initial value for the DDAE
[~,~,~,A2,f2,mu,~,~] = regularize_strange_ldae(E,A,@(t)B(t)*phi(t-tau)+f(t),t0,options);
x0 = x0 - pinv(A2)*(A2*x0+f2);

[t,xx] = solve_ldae_bvp(EE,AA,ff,eye(L*n),-diag(ones((L-1)*n,1),-n),[x0;zeros((L-1)*n,1)],[t0,t0+tau],options_new);

x = zeros(n,L*N+1);
x(:,1) = xx(1:n,1); 
for i = 1:L
    x(:,(i-1)*N+2:i*N+1) = xx((i-1)*n+1:i*n,2:N+1);
    t((i-1)*N+2:i*N+1) = t(2:N+1)+(i-1)*tau;
end

function EE = E_reformulated(E,m,n,L,t,tau)
EE = zeros(n*L);
for i=1:L
    EE((i-1)*m+1:i*m,(i-1)*n+1:i*n) = E(t+(i-1)*tau);
end

function AA = A_reformulated(A,B,m,n,L,t,tau)
AA = zeros(n*L);
AA(1:m,1:n) = A(t);
for i=2:L
    AA((i-1)*m+1:i*m,(i-1)*n+1:i*n) = A(t+(i-1)*tau);
    AA((i-1)*m+1:i*m,(i-2)*n+1:(i-1)*n) = B(t+(i-1)*tau);
end

function ff = f_reformulated(f,B,phi,m,n,L,t,tau)
ff = zeros(m*L,1);
ff(1:m) = f(t)+B(t)*phi(t-tau);
for i=2:L
    ff((i-1)*m+1:i*m) = f(t+(i-1)*tau);
end