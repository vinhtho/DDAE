function [E_1,A_1,f_1,A_2,f_2,mu,Z1,Z2] = regularize_strange_ldae(E,A,f,ti,options)
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

% set the options
isConst = 0;
KMax = 3;
mu0 = 0;
muMax = 3;
tolA = 1e-7;
tolR = 1e-7;
if exist('options','var')
    if isfield(options,'isConst')   isConst=options.isConst; end
    if isfield(options,'AbsTol')    tolA=options.AbsTol; end
    if isfield(options,'MaxStrIdx') muMax=options.MaxStrIdx; end
    if isfield(options,'RelTol')    tolR=options.RelTol; end
    if isfield(options,'StrIdx')    mu0=options.StrIdx; muMax = max(muMax,mu0); end
end
% tolerance for the matrix differential
tol=tolR;

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
    g = inflatef(f,ti,mu,tol,m);
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
%   END OF LINEAR_DDAE_INDEX_REDUCTION

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