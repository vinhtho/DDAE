function NM  =  inflateEA( E,A,t,mu,tolR )
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
end