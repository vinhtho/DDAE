%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   matrix_differential.m
%   ---------------------
%   
%   Approximates numerical the k-th-differential of a matrix-valued function A(s) w.r.
%   to s at the time point t.
%
%   INPUT PARAMETERS
%   ----------------
%   A   a matrix-valued function depending on a one dimensional parameter
%   t   the time for which the derivative is calculated
%   k   the order of the derivative
%   tol a measure for the quality of the approximation. The approximated differential 
%       will be returned, when the norm between two successive members of 
%       the approximation sequence is smaller than tol.
%
%   OUTPUT PARAMETERS
%   -----------------
%   dA  A approximation of the k-th derivative of A at the point t.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

end

