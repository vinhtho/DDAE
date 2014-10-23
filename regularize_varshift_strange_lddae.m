function [E_1,A_1,B_1,f_1,A_2,B_2,f_2,mu,K,Z1,Z2,U1] = regularize_varshift_strange_lddae(E,A,B,f,tau,ti,options)
%regularize_shift_strange_lddae
%
%   Subroutine for regularizing the DDAE
%           .
%       E(t)x(t) = A(t)x(t) + B(t)x(t-tau) + f(t),     t0 <= t <= tf,
%           x(t) = phi(t),                             t <= t0
%
%   with shift index bigger zero and strangeness index bigger zero.
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
%   B       fcn_handle      n-by-n matrix function
%   f       fcn_handle      n-by-1 inhomogeneity function
%   tau     double          the delay
%   xtau    fcn_handle      xtau(t) = x(t-tau)
%   ti      double          a time point
%   options struct          see other codes for explanation

% set the options
isConst = 0;
K0 = 0;
KMax = 3;
mu0 = 0;
muMax = 3;
tolA = 1e-7;
tolR = 1e-7;
if exist('options','var')
    if isfield(options,'isConst')   isConst=options.isConst; end
    if isfield(options,'AbsTol')    tolA=options.AbsTol; end
    if isfield(options,'MaxShift')  KMax=options.MaxShift; end
    if isfield(options,'MaxStrIdx') muMax=options.MaxStrIdx; end
    if isfield(options,'RelTol')    tolR=options.RelTol; end
    if isfield(options,'Shift')     K0=max(0,options.Shift); end
    if isfield(options,'StrIdx')    mu0=options.StrIdx; end
    if isfield(options,'x0')        x0=options.x0; end
end
% tolerance for the matrix differential
tol= tolR;
E0 = E(ti);
[m,n]=size(E0);
muMax = max(mu0,muMax);
%
% container for the derivative arrays of the original and timeshifted
% systems
NMP = zeros((muMax+1)*m*(KMax+1),(muMax+2)*n*(KMax+1));
%                                                      .
% logical indices corresponding to the coefficients of x and x in NMP
idx2M = false((muMax+2)*n*(KMax+1),1);
idx2M(n+1:2*n) = true;
idx2N = false((muMax+2)*n*(KMax+1),1);
idx2N(1:n) = true;
t_shifted = zeros(KMax+1,1);
t_shifted(1) = ti;
for i = 1:KMax
    t_shifted(i+1) = fsolve(@(t1) t1-tau(t1)-t_shifted(i),t_shifted(i)+tau(t_shifted(i)),optimset('Display','off'));
end
for K = 0:KMax
    % logical indices for selecting certain rows and columns in NMP
    idx1 = false((muMax+1)*m*(KMax+1),1);
    idx2 = false((muMax+2)*n*(KMax+1),1);
    
    for mu = 0:muMax
        NMP((1:(mu+1)*m)+K*(muMax+1)*m,(1:(mu+2)*n)+K*(muMax+2)*n) = inflateEA(E,A,t_shifted(K+1),mu,tolR);
        if K>0
            NMP((1:(mu+1)*m)+K*(muMax+1)*m,(1:(mu+1)*n)+(K-1)*(muMax+2)*n) = -inflateB_varshift(B,t_shifted(K+1),tau,mu,tolR,m,n);
            if K<K0
                continue;
            end
            % logical indices for the coefficients of x(t+tau), x(t+2*tau),
            % ..., x(t+K*tau) and its derivatives up to order mu in NMP
            for i = 1:K
                idx1((1:(mu+1)*m)+i*(muMax+1)*m)=true;
                idx2((1:(mu+2)*n)+i*(muMax+2)*n)=true;
            end
        end
        if mu<mu0
            continue;
        end
        % logical indices for the coefficients of x and its derivatives up 
        % to order mu in NMP
        idx1(1:(mu+1)*m) = true;
        idx2((2*n+1):(mu+2)*n) = true;
        %
        % extract a system without x(t+tau), x(t+2*tau), ..., x(t+K*tau) and its derivatives
        U1 = null2(NMP(idx1,idx2)',tolR);
        M = U1'*NMP(idx1,idx2M);
        N = -U1'*NMP(idx1,idx2N);
        P = U1(1:(mu+1)*m,:)'*inflateB_varshift(B,ti,tau,mu,tolR,m,n);
        %
        % TODO SOME COMMENTS
        Z2 = null2([M,P(:,n+1:(mu+1)*n)]',tolR);
        Z1 = orth2(M,tolR);
        %
        % extract the coefficients of the algebraic variables
        A_2 = Z2'*N;
        
        T2= null2(A_2,tolR);
        
        % check if the number of (linearly independent) algebraic equations
        % a and differential equations d is equal to the number of
        % variables n, if not then continue by increasing mu or K
        a = rank(A_2,tolR);
        if norm(isnan(Z1'*M*T2)+isinf(Z1'*M*T2),1)>0
            Z1(isnan(Z1))=0;
        end
        d = rank(Z1'*M*T2,tolR);
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
        B_2 = Z2'*P(:,1:n);
        
        % extract the coefficients of the differential variables
        E_1 = Z1'*M;
        % remove redundant equations
        Y1 = orth2(E_1*T2,tolR);
        E_1 = Y1'*E_1;
        % update the selector Z2
        Z1 = Z1*Y1;
        
        % check if the derivatives of x(t-tau) vanish in the differential
        % part
        B_1 = Z1'*P;
        if mu>0 && d>0
            if max(max(abs(B_1(:,n+1:end))))>tolR*max(max(B_1),1)
                mess1 = sprintf(['\nmaxmax(B_1(:,n+1:end))/max(norm(B_1,1),1)) = ',num2str(max(max(abs(B_1(:,n+1:end))))/max(norm(B_1,1),1))]);
                mess2 = sprintf(['\ntolerance                            = ',num2str(tolR)]);
                mess3 = sprintf('\n\nACCORDING TO THE CHOSEN TOLERANCE, THE SYSTEM IS OF ADVANCED TYPE, USING THE METHOD OF STEPS MIGHT PRODUCE LARGE ERRORS.');
                warning([mess1,mess2,mess3])
            end
        end
        B_1 = B_1(:,1:n);
        
        % extract the algebraic and differential parts for f and the
        % differential parts for E, A and B
        g = U1'*inflatef(f,t_shifted,K,mu,tol,m);
        f_2 = Z2'*g;
        A_1 = Z1'*N;
        f_1 = Z1'*g;
        return
    end
end

error('MAXIMAL NUMBER OF SHIFTS AND STRANGENESS REDUCTION STEPS REACHED. REGULARIZATION OF THE DDAE FAILED.')
%   END OF LINEAR_DDAE_INDEX_REDUCTION

function P = inflateB_varshift(B,tK,tau,mu,tol,m,n)
% builds the matrix
%    _                                          _
%   |                                            |
%   |   B                                        |
%   |   .          .                             |
%   |   B    B*(1-tau)                           |
%   |   ..    .    .     ..          .           |
%   |   B    2*B*(1-tau)-B*tau    B*(1-tau)^2    |
%   |                                            |
%   |_  etc.                                    _|
%
% i.e. the derivative array of B(t)*x(t-tau(t))
%
% it is hard coded up to mu = 3 because we haven't found a nice forumla yet
B0 = B(tK);
B1 = matrix_differential(B,tK,1,tol,m,n);
B2 = matrix_differential(B,tK,2,tol,m,n);
B3 = matrix_differential(B,tK,3,tol,m,n);
tau1 = matrix_differential(tau,tK,1,tol,1,1);
tau2 = matrix_differential(tau,tK,2,tol,1,1);
tau3 = matrix_differential(tau,tK,3,tol,1,1);
switch mu
    case 0
        P = B(tK);
    case 1
        P = [
            B0,zeros(m,n);
            B1,B0*(1-tau1)
            ];
    case 2
        P = [
            B0,zeros(m,2*n);
            B1,B0*(1-tau1),zeros(m,n);
            B2,2*B1*(1-tau1)-B0*tau2,B0*(1-tau1)^2
            ];
    case 3
        P = [
            B0,zeros(m,3*n);
            B1,B0*(1-tau1),zeros(m,2*n);
            B2,2*B1*(1-tau1)-B0*tau2,B0*(1-tau1)^2,zeros(m,n);
            B3,3*B2*(1-tau1)-3*B1*tau2-B0*tau3,3*B1*(1-tau1)^2-3*B0*(1-tau1)*tau2,B0*(1-tau1)^3
            ];
end


function g = inflatef(f,t_shifted,K,mu,tol,m)
%
% builds the vector
%    _                  _
%   |                    |
%   |   f (t_shifted(1)) |
%   |   .                |
%   |   f (t_shifted(1)) |
%   |   .                |
%   |   .                |
%   |   .                |
%   |    (mu)            |
%   |   f (t_shifted(1)) |
%   |   .                |
%   |   .                |
%   |   .                |
%   |    (mu)            |
%   |   f (t_shifted(K)) |
%   |_                  _|
%
g = zeros((K+1)*(mu+1)*m,1);
for j = 0:K
    for i = 1:(mu+1)
        g(((i-1)*m+1:i*m)+j*(mu+1)*m) = matrix_differential(f,t_shifted(j+1),i-1,tol,m,1);
    end
end