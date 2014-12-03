function [t,x,info] = solve_causal_ddae(E,A,B,f,tau,phi,tspan,options)
%SOLVE_CAUSAL_DDAE numerical solver for non-causal linear delay
% differential-algebraic equations of the form
%   E(t)\dot{x}(t)=A(t)x(t)+sum_i B_i(t)x(t-tau_i(t))+f(t)  for t\in(t0,tf]
%             x(t)=phi(t)                                   for t<=t0
% with multiple delay functions tau_i(t)>=tau_min>0 and history function 
% phi.
%
% The corresponding DAE Ex = Ax can have strangeness index bigger than
% zero (differentiation index bigger than one). However, note that bigger
% errors might occur if the strangeness index is too big, because
% derivatives are approximated by finite differences.
% The time stepping is based on the method of steps and using Radau IIA 
% collocation using constant step size.
%
% @parameters:
%   E,A,B       Coefficients of the DDAE, m-by-n matrix functions.
%   f           m-by-1 vector function.
%   tau         Variable lag, scalar function.
%   phi         n-by-1 history function.
%   tspan       Considered time interval [t0,tf].
%   options     Struct for optional parameters, set by
%               'options.FieldName = FieldValue', see below
%
% @options
%   Iter        Total number of time steps, default: 100.
%   Step        Step size.
%
%   AbsTol      Absolute tolerance, default: 1e-5.
%   RelTol      Relative tolerance, default: 1e-5.
%
%   StrIdx      Lower bound for the strangeness index.
%   MaxStrIdx   Upper bound for the strangeness index.
%
%   InitVal     Initial value, not necessarily consistent.
%   IsConst     Are E and A constant matrices?
%
% @supporting functions:
%   getRegularizedSystem
%   inflateEA
%   inflateB
%   inflatef
%   eval_all_hist_funcs
%   nevilleAitken
%
% @return values:
%   t           t(i) = t0+h_i with h_i the i-th step size.
%   x           numerical solution at the time nodes in t.
%   info        Struct with information.
%
% @author:
%       Vinh Tho Ma, TU Berlin, mavinh@math.tu-berlin.de
%       Phi Ha, TU Berlin, ha@math.tu-berlin.de

% set the options
%-------------------------------------------------------------------------%
% set missing fields in options
%-------------------------------------------------------------------------%
if ~exist('options','var'),options = {}; end

% number of iterations and step size (if both are set, step size has
% priority!)
if isfield(options,'Step')      
    options.Iter = diff(tspan)/options.Step; 
else
    options.Step = diff(tspan)/options.Iter;
end
if isfield(options,'Iter')  
    options.Step=diff(tspan)/options.Iter;
else
    options.Iter = 100;
end

% tolerances
if ~isfield(options,'AbsTol')    options.AbsTol = 1e-5; end
if ~isfield(options,'RelTol')    options.RelTol = 1e-5; end

% the DDAE's indeces (guesses, if unknown)
if ~isfield(options,'StrIdx')    options.StrIdx = 0; end
if ~isfield(options,'MaxStrIdx') options.MaxStrIdx = 3; end

% initial value (not necessarily consistent)
if ~isfield(options,'InitVal')   options.InitVal = phi(tspan(1)); end

% Are E and A constant matrices?
if ~isfield(options,'IsConst')   options.IsConst = false; end

%-------------------------------------------------------------------------%
% defining some parameters
%-------------------------------------------------------------------------%
t0 = tspan(1);
if not(isa(E,'function_handle'))
    error('E must be a function handle.'); 
end
[m,n] = size(E(0));
h = options.Step;
N = options.Iter;
tolR = options.RelTol;
x0 = options.InitVal;
isConst = options.IsConst;

% predefining info's fields
info.Solver = 'solve_causal_ddae';
info.Strangeness_index = -1;
info.Advanced = 0;
info.Number_of_differential_eqs = -1;
info.Number_of_algebraic_eqs = -1;
info.Computation_time = -1;

%-------------------------------------------------------------------------%
% some more input checks
%-------------------------------------------------------------------------%
% checking tau
if isa(tau,'double')
    tau = @(t)tau;
end
if not(isa(tau,'function_handle'))
    error('Delay tau must be a function handle.'); 
end
l = numel(tau(0));
% checking A
if not(isa(A,'function_handle'))
    error('A must be a function handle.'); 
end
if or(size(A(0),1)~=m,size(A(0),2)~=n)
    error('A has wrong size.')
end
% checking B
if not(isa(B,'function_handle'))
    error('B must be a function handle.'); 
end
if or(size(B(0),1)~=m,size(B(0),2)~=l*n)
    error('B has wrong size.')
end
% checking f
if not(isa(f,'function_handle'))
    error('f must be a function handle.'); 
end
if or(size(f(0),1)~=m,size(f(0),2)~=1)
    error('f has wrong size.')
end
% checking phi
if or(size(phi(0),1)~=n,size(phi(0),2)~=1)
    error('phi has wrong size.')
end
if options.MaxStrIdx<options.StrIdx
    error('MaxStrIdx must not be less than StrIdx.')
end

% the data for the RADAU IIA collocation, V is the inverse of A in the
% Butcher tableau
c=[(4-sqrt(6))/10; (4+sqrt(6))/10; 1];
V=[ 3.224744871391589   1.167840084690405  -0.253197264742181
    -3.567840084690405   0.775255128608412   1.053197264742181
    5.531972647421811  -7.531972647421810   5.000000000000000 ];
v0=-V*ones(3,1);

% the container for the approximate solution of the DDAE and its stage values
x=nan(3*n,N+1);
x(1:n,1)=phi(t0+(c(1)-1)*h);
x(n+1:2*n,1)=phi(t0+(c(2)-1)*h);

% find nearest consistent initial value to x0 by determining local
% strangeness-free form at t=t0 and replace x0 by it
hist_funcs = cell(l,1);
for k=1:l
    hist_funcs{k} = phi;
end
BXTAUF = @(s) B(s)*eval_all_hist_funcs(hist_funcs,s,tau(s),n,l)+f(s);
[E1,A1,~,A2,g2,mu,Z1,Z2] = getRegularizedSystem(E,A,BXTAUF,t0,options);
info.Strangeness_index = mu;
info.Number_of_differential_eqs = size(E1,1);
info.Number_of_algebraic_eqs = size(A2,1);

% compute nearest consistent initial vector to x0
x(2*n+1:3*n,1)=x0-pinv(A2)*(A2*x0+g2);

% the discretized time interval
t=t0+h*(0:1:N);

% containers for the matrices of the local strangeness-free formulations at
% t_ij = T(i)+h*c(j)
Etij=nan(n,n,3);
Atij=nan(n,n,3);
ftij=nan(n,1,3);

%-------------------------------------------------------------------------%
% Time integration.
%-------------------------------------------------------------------------%
tic
for i=1:N
    % For each collocation point...
    for j=1:3
        TAU = tau(t(i)+c(j)*h);
        hist_funcs = cell(1,l);
        for k = 1:l
            % determine "x(t-tau)" by using the histpry function phi or by
            % using interpolation
            if t(i)+c(j)*h-TAU(k)<t0
                hist_funcs{k} = phi;
            else
                % find the biggest time node smaller than t_i+c_j*h-tau
                t_tau_index=find(t(i)+c(j)*h-TAU(k)<t(1:i),1)-1;
                if isempty(t_tau_index)
                    warning('ONE DELAY BECAME SMALLER THAN THE STEP SIZE. LONG STEPS NOT IMPLEMENTED YET IN solve_''causal_ddae.m''. TERMINATING SOLVING PROCESS.')
                    t=t(1:i);
                    x=x(2*n+1:3*n,1:i);
                    info.Computation_time = toc;
                    return;
                end 
                t_tau=t(t_tau_index);
                % if t_i+c_j*h-tau is not a node point, i.e. not in t, then we
                % have to interpolate
                % prepare some data for the interpolation
                % we use a polynomial of degree 3, so we need 4 data points
                x0_tau=x(2*n+1:3*n,t_tau_index);
                X_tau=reshape(x(:,t_tau_index+1),n,3);
                hist_funcs{k} = @(s) nevilleAitken(t_tau+[0;c]*h,[x0_tau,X_tau],s);
            end
        end
        
        BXTAUF = @(s) B(s)*eval_all_hist_funcs(hist_funcs,s,tau(s),n,l)+f(s);
        
        % calculate locally regularized (i.e. strangeness-free) form at t =
        % t(i)-c(j)*h
        % we already have E1 and A1, if isConst is TRUE
        if isConst
            g = zeros((mu+1)*m,1);
            for k = 0:mu
                g(k*m+1:(k+1)*m) = matrix_differential(BXTAUF,t(i)+c(j)*h,k,tolR,m,1);
            end
            g1 = Z1'*(BXTAUF(t(i)+c(j)*h));
            if numel(Z2)>0
                g2 = Z2'*g;
            else
                g2 = zeros(0,1);
            end
        else
            [E1,A1,g1,A2,g2] = getRegularizedSystem(E,A,BXTAUF,t(i)+c(j)*h,options);
        end
        
        % check if the derivatives of x(t-tau) vanish in the differential
        % part, if not already known
        if info.Advanced == false
            P = inflateB(B,t(i)+c(j)*h,mu,tolR,m,l*n);
            if numel(Z2)>0
                B_2 = Z2'*P;
                if mu>0
                    if (max(max(abs(B_2(:,l*n+1:end))))>tolR*max(max(B_2),1))
                        warning('ACCORDING TO THE CHOSEN TOLERANCE, THE SYSTEM IS VERY LIKELEY OF ADVANCED TYPE, USING THE METHOD OF STEPS MIGHT PRODUCE LARGE ERRORS.')
                        info.Advanced = true;
                    end
                end
            end
        end
        Etij(:,:,j)=[E1;zeros(size(A2))];
        Atij(:,:,j)=[A1;A2];
        ftij(:,j)=[g1;g2];
    end
    
    % SOLVING THE LINEAR SYSTEM
    AA=zeros(3*n);
    bb=zeros(3*n,1);
    
    for j=1:3
        for k=1:3
            AA((j-1)*n+1:j*n,(k-1)*n+1:k*n)=Etij(:,:,j)/h*V(j,k)-(j==k)*Atij(:,:,j);
        end
        bb((j-1)*n+1:j*n)=ftij(:,j)-Etij(:,:,j)/h*v0(j)*x(2*n+1:3*n,i);
    end
    
    % the solution is a vector with length 3*n, it consists of the 3 values
    % of the polynomial at the collocation points T(i)+c(j)*h, j=1..3
    x(:,i+1)=AA\bb;
    
end


% "cutting out" the approximate solution
x=x(2*n+1:3*n,:);
info.Computation_time = toc;

function [E_1,A_1,f_1,A_2,f_2,mu,Z1,Z2] = getRegularizedSystem(E,A,f,ti,options)
%regularize_strange_ldae
%
%   Subroutine for regularzing the DDAE
%           .
%       E(t)x(t) = A(t)x(t) + B(t)x(t-tau) + f(t),     t0 <= t <= tf,
%           x(t) = phi(t),                             t <= t0.
%
%   The index reduction procedure was taken from
%   -----------------------------------------------------------------------
%   P. Kunkel, V. Mehrmann: Differential-Algebraic Equations, Chapter 6.1.
%   -----------------------------------------------------------------------
%
%   INPUT   DATA TYPE       COMMENTS
%   -----   ---------       --------
%   E       fcn_handle      n-by-n leading matrix function
%   A       fcn_handle      n-by-n matrix function
%   f       fcn_handle      n-by-1 inhomogeneity function
%   ti      double          a time point
%   options struct          see other codes for explanation
%

 muMax=options.MaxStrIdx; 
 tolR=options.RelTol;  
 mu0=options.StrIdx; 
 muMax = max(muMax,mu0);

E0 = E(ti);
[m,n]=size(E0);
muMax = max(mu0,muMax);

for mu = mu0:muMax
    NM = inflateEA(E,A,ti,mu,tolR);
    M = NM(:,(n+1):end);
    N = -NM(:,1:n);
    
    % TODO SOME COMMENTS
    Z2 = null2(M',tolR);
    %
    % extract the coefficients of the algebraic variables
    A_2 = Z2'*N;
    
    T2 = null2(A_2,tolR);
    
    % check if the number of (linearly independent) algebraic equations
    % a and differential equations d is equal to the number of
    % variables n, if not then continue by increasing mu or K
    a = rank(A_2,tolR);
    
    d = rank(E0*T2,tolR);
    if a+d~=n
        continue
    end
    
    % remove redundant algebraic equations
    if size(A_2,1)>0
        Y2 = orth2(A_2,tolR);
        A_2 = Y2'*A_2;
        % update the selector Z2
        Z2 = Z2*Y2;
    end
    
    % remove redundant equations
    Z1 = orth2(E0*T2,tolR);
    E_1 = Z1'*E0;
    
    % extract the algebraic and differential parts for f and the
    % differential parts for E, A and B
    g = inflatef(f,ti,mu,tolR,m);
    if a>0
        f_2 = Z2'*g;
    else
        f_2 = zeros(0,1);
    end
    A_1 = Z1'*A(ti);
    f_1 = Z1'*f(ti);
    return
end
warning('MAXIMAL NUMBER OF SHIFTS AND STRANGENESS REDUCTION STEPS REACHED. REGULARIZATION OF THE DDAE FAILED.')

%-------------------------------------------------------------------------%
% Auxiliary functions used during the time integration
%-------------------------------------------------------------------------%
function XTAU = eval_all_hist_funcs(hist_funcs,s,TAU,n,l)
XTAU = zeros(l*n,1);
for k = 1:l
    XTAU((k-1)*n+1:k*n,1) = feval(hist_funcs{k},s-TAU(k));
end
function px = nevilleAitken(X,F,x)
n=length(X);
%at first px is a container for the values in the Newton scheme
px=F;
% beginning the Newton scheme, see Numerische Mathematik 1
for i=1:n-1
    for j=1:n-i
        px(:,j)=((x-X(j))*px(:,j+1)-(x-X(j+i))*px(:,j))/(X(j+i)-X(j));
    end
end
px=px(:,1);
%-------------------------------------------------------------------------%
% Auxiliary functions for getRegulizedSystem()
%-------------------------------------------------------------------------%
function NM = inflateEA( E,A,t,mu,tolR )
%   INFLATEEA: Computes the derivative array of (E,A).
%
%   This file computes the derivative array of the homogeneous DAE
%                           .
%                       E(t)x(t) = A(t)x(t)
%
%   by differentiating mu times.
%
%   INPUT   
%   -----   
%   E       fcn_handle      m-by-n leading matrix function
%   A       fcn_handle      m-by-n matrix function
%   t       double          the time
%   mu      double          the strangeness index
%   tolR    double          the relative tolerance
%
%
%   OUTPUT      
%   ------
%   NM = [-N,M], where
%       M       double((mu+1)*m,(mu+1)*n)   M like in [1]
%       N       double((mu+1)*m,n)          first n columns of N in [1]
%
%   References:
%       [1] P.Kunkel, V. Mehrmann: Differential-Algebraic Equations,
%           chapters 3.1 and 3.2
%

E0 = E(t);
A0 = A(t);

[m,n] = size(E0);

dE = zeros((mu+1)*m,n);
dA = zeros((mu+1)*m,n);

NM = zeros((mu+1)*m,(mu+2)*n);

dE(1:m,1:n) = E0;
dA(1:m,1:n) = A0;
NM(1:m,n+1:2*n) = E0;

for l = 1:mu
    % make dE and dA contain all derivatives up to order l
    dE(l*m+1:(l+1)*m,1:n) = matrix_differential( E,t,l,tolR,m,n);
    dA(l*m+1:(l+1)*m,1:n) = matrix_differential( A,t,l,tolR,m,n);
    
    %Expand M_(l-1) to M_l
    for j = 0:l-1
        k = l-j;
        NM(l*m+1:(l+1)*m,(j+1)*n+1:(j+2)*n) = nchoosek(l,j)*dE(k*m+1:(k+1)*m,:)-nchoosek(l,j+1) * dA((k-1)*m+1:k*m,:);
    end
    NM(l*m+1:(l+1)*m,(l+1)*n+1:(l+2)*n) = dE(1:m,1:n);
end
NM(:,1:n) = -dA;
function g = inflatef(f,ti,mu,tol,m)
%
% builds the vector
%    _        _
%   |          |
%   |   f(ti)  |
%   |   .      |
%   |   f(ti)  |
%   |   ..     |
%   |   f(ti)  |
%   |  ...     |
%   |   f(ti)  |
%   |   .      |
%   |   .      |
%   |   .      |
%   |    (mu)  |
%   |   f (ti) |
%   |_        _|
%
g = zeros((mu+1)*m,1);
    for i = 1:(mu+1)
        g(((i-1)*m+1:((i-1)+1)*m)) = matrix_differential(f,ti,i-1,tol,m,1);
    end
function P = inflateB(B,ti,mu,tol,m,n)
%
% builds the matrix
%    _                                         _
%   |                                           |
%   |   B                                       |
%   |   .                                       |
%   |   B    B                                  |
%   |   ..   .                                  |
%   |   B   2B    B                             |
%   |  ...   ..   .                             |
%   |   B   3B   3B    B                        |
%   |   .    .    .    .    .                   |
%   |   .    .    .    .    .    .              |
%   |   .    .    .    .    .    .    .         |
%   |    (mu)                                   |
%   |   B    .    .    .    .    .    .    B    |
%   |_                                         _|
%
P = zeros((mu+1)*m,(mu+1)*n);
for i = 1:(mu+1)
    P((i-1)*m+1:((i-1)+1)*m,1:n) = matrix_differential( B,ti,i-1,tol,m,n);
    for j=1:i-1
        k = i-j-1;
        P((i-1)*m+1:((i-1)+1)*m,j*n+1:(j+1)*n) = round(prod((i-j:i-1)./(1:j)))*P(k*m+1:(k+1)*m,1:n);
    end
end
function dA = matrix_differential(A,t,k,tol,m,n)
%Parameters 
eps=0.01;
j=0;
delta=sqrt(eps*max(0.01,abs(t)));

% [m,n]=size(A(t));
temp=zeros(m,n,k+1);
% dA=A(0);
alpha=tol+1;

while j<2 && alpha>tol
    delta=delta/2;
    dA_old=A(0);
    for i=0:k
%         temp(:,:,i+1)=(-1)^i*nchoosek(k,i)*A(t+(k/2-i)*delta);
        temp(:,:,i+1)=(-1)^i*round(prod(((k-i+1):k)./(1:i)))*A(t+(k/2-i)*delta);
    end
    dA=sum(temp,3)/delta^k;
    alpha=norm(dA-dA_old);
    j=j+1;
end

if min(min(isfinite(dA)))==0
    warning('ERROR IN maxtrix_differential.m!')
end
function Z = null2(A,tol)
[m,n] = size(A);
[U,S,V]=svd(A,0);
if m > 1
    s = diag(S);
  elseif m == 1
      s = S(1);
  else s = 0;
end
r = sum(s > max(m,n) * max(s(1),1) * tol);
Z = V(:,r+1:n);
function Q = orth2(A,tol)
if isempty(A)
    Q=A;
    return;
end
[U,S] = svd(A,0);
[m,n] = size(A);
if m > 1, s = diag(S);
   elseif m == 1, s = S(1);
   else s = 0;
end
r = sum(s > max(m,n) * max(s(1),1) * tol);
Q = U(:,1:r);