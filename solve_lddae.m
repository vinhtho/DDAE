function [t,x] = solve_lddae(E,A,B,f,tau,phi,tspan,options)
%solve_lddae
%
%   Solver for linear delay differential algebraic equations (DDAEs) with l
%   time-dependent delays tau_1(t),...,tau_l(t) > 0, i.e.
%       .
%   E(t)x(t) = A(t)x(t) + B(t)*[x(t-tau_1(t));...;x(t-tau_l(t))] + f(t),
%
%   for t0 <= t <= tf and x(t) = phi(t) for t <= t0,
%
%   based on using the method of steps and using Radau IIA collocation. The
%   delayed part and the inhomogeinity f are considered as an inhomgeinity
%                                       .
%   g and we solve the solve the "DAE" Ex = Ax + g step by step. This DAE
%   can have strangeness index bigger than zero (differentiation index
%   bigger than one). However, note that bigger errors might occur if the
%   strangeness index is too big, because derivatives are approximated by
%   finite differences.
%
%   INPUT
%   -----
%   E       fcn_handle      m-by-n matrix function
%   A       fcn_handle      m-by-n matrix function
%   B       fcn_handle      m-by-n-by-l matrix function
%   f       fcn_handle      m-by-1 inhomogeneity function
%   tau     fcn_handle      1-by-l delay function
%   phi     fcn_handle      n-by-1 history function
%   tspan   double(1,2)     the solution interval
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

%tau has to be a function, so we turn constant delays into a function
if not(isa(tau,'function_handle'))
    tau=@(s)tau;
end

% the initial time
t0 = tspan(1);

% numbers of delays
l=length(tau(t0));

E0 = E(tspan(1));

% the dimension of the state variable
[m,n] = size(E0);


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
[E1,A1,~,A2,g2,mu,Z1,Z2] = regularize_strange_ldae(E,A,BXTAUF,t0,options);

x(2*n+1:3*n,1)=x0-pinv(A2)*(A2*x0+g2);

% the discretized time interval
t=t0+h*(0:1:N);

% containers for the matrices of the local strangeness-free formulations at
% t_ij = T(i)+h*c(j)
Etij=nan(n,n,3);
Atij=nan(n,n,3);
ftij=nan(n,1,3);

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
                t_tau_index=find(t(i)+c(j)*h-TAU(k)<t,1)-1;
                t_tau=t(t_tau_index);
                % if t_i+c_j*h-tau is not a node point, i.e. not in t, then we
                % have to interpolate
                % prepare some data for the interpolation
                % we use a polynomial of degree 3, so we need 4 data points
                x0_tau=x(2*n+1:3*n,t_tau_index);
                X_tau=reshape(x(:,t_tau_index+1),n,3);
                hist_funcs{k} = @(s) poleval_neville_aitken(t_tau+[0;c]*h,[x0_tau,X_tau],s);
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
            g2 = Z2'*g;
        else
            [E1,A1,g1,A2,g2] = regularize_strange_ldae(E,A,BXTAUF,t(i)+c(j)*h,options);
        end
        
        P = inflateB(B,t(i)+c(j)*h,mu,tolR,m,l*n);
        % check if the derivatives of x(t-tau) vanish in the differential
        % part
        B_2 = Z2'*P;
        if mu>0
            if max(max(abs(B_2(:,l*n+1:end))))>tolR*max(max(B_2),1)
                mess1 = sprintf(['\nmaxmax(B_2(:,l*n+1:end))/max(norm(B_2,1),1)) = ',num2str(max(max(abs(B_2(:,n+1:end))))/max(norm(B_2,1),1))]);
                mess2 = sprintf(['\ntolerance                            = ',num2str(tolR)]);
                mess3 = sprintf('\n\nACCORDING TO THE CHOSEN TOLERANCE, THE SYSTEM IS OF ADVANCED TYPE, USING THE METHOD OF STEPS MIGHT PRODUCE LARGE ERRORS.');
                warning([mess1,mess2,mess3])
            end
        end
        
        Etij(:,:,j)=[E1;zeros(size(A2))];
        Atij(:,:,j)=[A1;A2];
        ftij(:,j)=[g1;g2];
    end
    
    % SOLVING THE LINEAR SYSTEM
    
    % Please have a look at page 4 of the paper above. There are 4
    % conditions a)-d). What is done here is the solving the system in
    % condition a), which is linear in our case.
    
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

% end of LDDAE_SOLVE

function XTAU = eval_all_hist_funcs(hist_funcs,s,TAU,n,l)
XTAU = zeros(l*n,1);
for k = 1:l
    XTAU((k-1)*n+1:k*n,1) = feval(hist_funcs{k},s-TAU(k));
end

function px = poleval_neville_aitken(X,F,x)

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