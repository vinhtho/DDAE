function [t,x] = solve_ldae_bvp(E,A,f,C,D,r,tspan,options)
%solve_ldae_bvp
%
%   Solver for a linear DAE
%       .
%   E(t)x(t) = A(t)x(t) + f(t)      for t in tspan = [t0;tN]
%
%   with linear boundary conditions
%
%   C*x(t0) + D*x(tf) = r.
%
%   We use the Radau5 collocation method for the discretization.
%   
%   The strangeness index of the DAE can be bigger than zero 
%   (differentiation index bigger than one).
%
%   INPUT   
%   -----   
%   E       fcn_handle      m-by-n leading matrix function
%   A       fcn_handle      m-by-n matrix function
%   f       fcn_handle      m-by-1 inhomogeneity function
%   C,D     double(>d,n)    d is the number of differential equations in
%                           the strangeness-free reformulation of the DAE
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
if exist('options','var')
    if isfield(options,'isConst')   isConst=options.isConst; end
    if isfield(options,'AbsTol')    tolA=options.AbsTol; end
    if isfield(options,'MaxStrIdx') muMax=options.MaxStrIdx; end
    if isfield(options,'NGrid')     N=options.NGrid-1; h=diff(tspan)/N; end
    if isfield(options,'RelTol')    tolR=options.RelTol; end
    if isfield(options,'StepSize')  h=options.StepSize; N=floor(diff(tspan)/h); end
    if isfield(options,'StrIdx')    mu0=options.StrIdx; muMax = max(muMax,mu0); end
else
    options = {};
end


c=[(4-sqrt(6))/10; (4+sqrt(6))/10; 1];

V=[ 3.224744871391589   1.167840084690405  -0.253197264742181
   -3.567840084690405   0.775255128608412   1.053197264742181
    5.531972647421811  -7.531972647421810   5.000000000000000 ];

v0=-V*ones(3,1);

% the dimension
[m,n] = size(E(tspan(1)));

% the stepsize
h = diff(tspan)/N;

% the discretized time
t=tspan(1):h:tspan(2);

% containers for the big linear system
AA = sparse(zeros((3*N+1)*n));
bb = zeros((3*N+1)*n,1);

for i = 1:N
    %disp(['Reducing index at t = ',num2str(t(i)+c(1)*h),', ',num2str(t(i)+c(2)*h),' and ',num2str(t(i)+c(3)*h)])
    for j = 1:3
        
        if isConst && exist('Z1','var')
            f1 = Z1'*f(t(i)+c(j)*h);
            f2 = Z2'*inflatef(f,t(i)+c(j)*h,S,tolR,m);
        else
            [E1,A1,f1,A2,f2,S,Z1,Z2] = regularize_strange_ldae(E,A,f,t(i)+c(j)*h,options);
        disp(['Strangeness-index S = ',num2str(S)])
            
        end
        
        
        Eij = [E1;zeros(size(A2))];
        Aij = [A1;A2];
        fij = [f1;f2];
        
        AA((3*(i-1)+j-1)*n+1:(3*(i-1)+j)*n,3*(i-1)*n+1:(3*(i-1)+1)*n) = Eij*v0(j)/h;
        
        for l = 1:3
            AA((3*(i-1)+j-1)*n+1:(3*(i-1)+j)*n,(3*(i-1)+l)*n+1:(3*(i-1)+l+1)*n) = Eij*V(j,l)/h - (l==j)*Aij;
        end
        
        bb((3*(i-1)+j-1)*n+1:(3*(i-1)+j)*n) = fij;
        
    end
end

% the boundary conditions
AA(3*N*n+1:(3*N+1)*n,1:n) = C;
AA(3*N*n+1:(3*N+1)*n,3*N*n+1:(3*N+1)*n) = D;
bb(3*N*n+1:(3*N+1)*n) = r;

disp(['Solving linear system of dimension n=',num2str((3*N+1)*n)])

xx = AA\bb;

% reshaping the solution to an appropriate size
x = reshape(xx(n+1:end),3*n,N);
x = [xx(1:n),x(2*n+1:3*n,:)];

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