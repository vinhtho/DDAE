function [t,x] = solve_varshifted_lddae(E,A,B,f,tau,phi,tspan,options)
%solve_varshifted_lddae
%
%   Solver for linear delay differential algebraic equations (DDAEs) with
%   variable coefficients and single, variable delay tau > 0, i.e.
%           .
%       E(t)x(t) = A(t)x(t) + B(t)x(t-tau) + f(t),     t0 <= t <= tf,
%           x(t) = phi(t),                             t <= t0,
%
%   based on using the method of steps and using Radau IIA collocation.
%                          .
%   The corresponding DAE Ex = Ax can have strangeness index bigger than
%   zero (differentiation index bigger than one). However, note that bigger
%   errors might occur if the strangeness index is too big, because
%   derivatives are approximated by finite differences.
%
%   Furthermore, we allow the so called shift index to be bigger than zero.
%
%   A regularization routine regularizes the DDAE by shifting appropriate
%   equations and reducing the strangeness.
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
%   MaxShift    double      upper bound for the shift index
%   MaxStrIdx   double      upper bound for the strangeness index
%   NGrid       double      no. of grid points for time discretization
%   RelTol      double      relative tolerance
%   Shift       double      shift index (if known)
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

%-------------------------------------------------------------------------%
% set the options
%-------------------------------------------------------------------------%
isConst = 0;
K0 = 0;
KMax = 3;
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
    if isfield(options,'MaxIter')   MaxIter = options.MaxIter; end
    if isfield(options,'MaxShift')  KMax=options.MaxShift; end
    if isfield(options,'MaxStrIdx') muMax=options.MaxStrIdx; end
    if isfield(options,'NGrid')     N=options.NGrid-1; h=diff(tspan)/N; end
    if isfield(options,'RelTol')    tolR=options.RelTol; end
    if isfield(options,'Shift')     K0=max(1,options.Shift); end
    if isfield(options,'StepSize')  h=options.StepSize; N=floor(diff(tspan)/h); end
    if isfield(options,'StrIdx')    mu0=options.StrIdx; end
    if isfield(options,'x0')        x0=options.x0; end
else
    options = {};
end

%-------------------------------------------------------------------------%
% doing something with the input
%-------------------------------------------------------------------------%
t0 = tspan(1);
E0 = E(tspan(1));
[m,n] = size(E0);


%-------------------------------------------------------------------------%
% preparation for the main loop
%-------------------------------------------------------------------------%
% The vector c comes from the Butcher tableau of the 3-stage Radar IIa
% method.
c=[(4-sqrt(6))/10; (4+sqrt(6))/10; 1];

% The discretized time interval (equidistant at the moment)
t=t0+h*(0:1:N);

% The column vector x(:,j) with length 3*n contains an approximation of the
% DDAE's solution at the time points t(j-1)+c(1)*h, t(j-1)+c(2)*h and t(j)
% in the first n, second n and third n entries, respectively, where h
% denotes the (current) step size.
x=nan(3*n,N+1);
x(1:n,1)=phi(t0+(c(1)-1)*h);
x(n+1:2*n,1)=phi(t0+(c(2)-1)*h);


%-------------------------------------------------------------------------%
% finding consistent initial value
%-------------------------------------------------------------------------%
[E1,A1,B1,~,A2,B2,f2,mu,K,Z1,Z2,U1] = regularize_varshift_strange_lddae(E,A,B,f,tau,t0,options);
x(2*n+1:3*n,1)=x0-pinv(A2)*(A2*x0+B2*phi(t0-tau(t0))+f2);

%-------------------------------------------------------------------------%
% main loop - time integration
%-------------------------------------------------------------------------%
for i=1:N
    % Compute approximation of x at t = t(i)+h.
    x(:,i+1) = performLongStep(E,A,B,f,tau,phi,t(1:i),x(:,1:i),h,options);
    t(i+1) = t(i) + h;
end

% "cutting out" the approximate solution at t
x=x(2*n+1:3*n,:);

%-------------------------------------------------------------------------%
% perform one step of the method of steps, h is smaller than tau
%-------------------------------------------------------------------------%
function x_next = performNormalStep(E,A,B,f,tau,phi,t,x,h,options)
if h>tau(t(end)+h)
    error('Step size is smaller then the delay.')
end
n=numel(x(:,1))/3;
% The vector c comes from the Butcher tableau of the 3-stage Radar IIa
% method.
c=[(4-sqrt(6))/10; (4+sqrt(6))/10; 1];
% containers for the matrices of the local strangeness-free formulations at
% t_ij = T(i)+h*c(j)
Etij=nan(n,n,3);
Atij=nan(n,n,3);
ftij=nan(n,3);
xtau=nan(n,3);
% For each collocation point...
for j=1:3
    % determine "x(t-tau)"
    t_tau = t(end)+c(j)*h-tau(t(end)+c(j)*h);
    if t_tau<t(1)
        % t-tau is less than t0, use the history function
        xtau(:,j) = phi(t_tau);
    else
        % t-tau is in the interval [t(1),t(end)]
        % find the biggest time node smaller than t_i+c_j*h-tau
        L = find(t_tau<t,1)-1;
        % if t_i+c_j*h-tau is not a node point, i.e. not in t, then we
        % have to interpolate
        % prepare some data for the interpolation
        % we use a polynomial of degree 3, so we need 4 data points
        x0_tau=x(2*n+1:3*n,L);
        X_tau=reshape(x(:,L+1),n,3);
        h_tau = t(L+1)-t(L);
        % interpolate with Neville-Aitken
        xtau(:,j) = poleval_neville_aitken(t(L)+[0;c]*h_tau,[x0_tau,X_tau],t_tau);
    end
    % calculate locally regularized form at t = t(i)-c(j)*h
    [E1,A1,B1,f1,A2,B2,f2] = regularize_varshift_strange_lddae(E,A,B,f,tau,t(end)+c(j)*h,options);
    Etij(:,:,j)=[E1;zeros(size(A2))];
    Atij(:,:,j)=[A1;A2];
    ftij(:,j)=[B1*xtau(:,j)+f1;B2*xtau(:,j)+f2];
end
% SOLVING THE LINEAR SYSTEM
x_next = solveLinSystem(Etij,Atij,ftij,x,h);


%-------------------------------------------------------------------------%
% perform a long step, i.e. h>tau and x(t-tau) has to be predicted using
% extrapolation of the last computed cubic polynomial
%-------------------------------------------------------------------------%
function x_next = performLongStep(E,A,B,f,tau,phi,t,x,h,options)
n=numel(x(:,1))/3;
% The vector c comes from the Butcher tableau of the 3-stage Radar IIa
% method.
c=[(4-sqrt(6))/10; (4+sqrt(6))/10; 1];
% containers for the matrices of the local strangeness-free formulations at
% t_ij = T(i)+h*c(j)
Etij=nan(n,n,3);
Atij=nan(n,n,3);
ftij=nan(n,3);
xtau=nan(n,3);
for i=1:10
    % For each collocation point...
    if i==1
    for j=1:3
        % determine "x(t-tau)"
        t_tau = t(end)+c(j)*h-tau(t(end)+c(j)*h);
        if t_tau<t(1)
            % t-tau is less than t0, use the history function
            xtau(:,j) = phi(t_tau);
        elseif t_tau<=t(end)
            % t-tau is in the interval [t(1),t(end)]
            % find the biggest time node smaller than t_i+c_j*h-tau
            L = find(t_tau<t,1)-1;
            % if t_i+c_j*h-tau is not a node point, i.e. not in t, then we
            % have to interpolate
            % prepare some data for the interpolation
            % we use a polynomial of degree 3, so we need 4 data points
            x0_tau=x(2*n+1:3*n,L);
            X_tau=reshape(x(:,L+1),n,3);
            h_tau = t(L+1)-t(L);
            % interpolate with Neville-Aitken
            xtau(:,j) = poleval_neville_aitken(t(L)+[0;c]*h_tau,[x0_tau,X_tau],t_tau);
        else
            x0_tau=x(2*n+1:3*n,end-1);
            X_tau=reshape(x(:,end),n,3);
            h_tau = t(end)-t(end-1);
            % interpolate with Neville-Aitken
            xtau(:,j) = poleval_neville_aitken(t(end-1)+[0;c]*h_tau,[x0_tau,X_tau],t_tau);
        end
    end
    end
    for j=1:3
        % calculate locally regularized form at t = t(i)-c(j)*h
        [E1,A1,B1,f1,A2,B2,f2] = regularize_varshift_strange_lddae(E,A,B,f,tau,t(end)+c(j)*h,options);
        Etij(:,:,j)=[E1;zeros(size(A2))];
        Atij(:,:,j)=[A1;A2];
        ftij(:,j)=[B1*xtau(:,j)+f1;B2*xtau(:,j)+f2];
    end
    % SOLVING THE LINEAR SYSTEM
    x_next = solveLinSystem(Etij,Atij,ftij,x,h);
    
    % correcting xtau
    xtau_new = xtau;
    for j=1:3
        t_tau = t(end)+c(j)*h-tau(t(end)+c(j)*h);
        x0_tau=x(2*n+1:3*n,end);
        X_tau=reshape(x_next,n,3);
        % interpolate with Neville-Aitken
        xtau_new(:,j) = poleval_neville_aitken(t(end)+[0;c]*h,[x0_tau,X_tau],t_tau);
    end
    
    if norm(xtau_new-xtau)<1e-5
        i
        return
    else
        xtau=xtau_new;
    end
end

%-------------------------------------------------------------------------%
% auxiliary function for solving the big linear system
%-------------------------------------------------------------------------%
function x_next = solveLinSystem(Etij,Atij,ftij,x,h)
n=size(Etij,1);
% The matrix V is the inverse of A in the Butcher tableau of the Radar IIa 
% method. 
V=[ 3.224744871391589   1.167840084690405  -0.253197264742181
    -3.567840084690405   0.775255128608412   1.053197264742181
    5.531972647421811  -7.531972647421810   5.000000000000000 ];
v0=-V*ones(3,1);
AA=zeros(3*n);
bb=zeros(3*n,1);
for j=1:3
    for k=1:3
        AA((j-1)*n+1:j*n,(k-1)*n+1:k*n)=Etij(:,:,j)/h*V(j,k)-(j==k)*Atij(:,:,j);
    end
    bb((j-1)*n+1:j*n)=ftij(:,j)-Etij(:,:,j)/h*v0(j)*x(2*n+1:3*n,end);
end
% the solution is a vector with length 3*n, it consists of the 3 values
% of the polynomial at the collocation points t(i)+c(j)*h, j=1..3
x_next=AA\bb;
 
%-------------------------------------------------------------------------%
% auxiliary function for interpolating/extrapolating
%-------------------------------------------------------------------------%
function px = poleval_neville_aitken(X,F,x)

n=length(X);
% at first px is a container for the values in the Newton scheme
px=F;
% beginning the Newton scheme, see Numerische Mathematik 1
for i=1:n-1
    for j=1:n-i
        px(:,j)=((x-X(j))*px(:,j+1)-(x-X(j+i))*px(:,j))/(X(j+i)-X(j));
    end
end
px=px(:,1);