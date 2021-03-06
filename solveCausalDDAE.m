function [t,x,info] = coldedaNonCausal(E,A,B,f,tau,phi,tspan,options)
%COLDEDA_NONCAUSAL numerical solver for non-causal linear delay
% differential-algebraic equations of the form
%   E(t)\dot{x}(t) = A(t)x(t) + B(t)x(t-tau(t)) + f(t)  for t\in(t0,tf]
%             x(t) = phi(t),                            for t<=t0
% with smooth delay tau(t)>=0 and history function phi.
%
% Noncausal means that the so called shift index can be bigger than zero,
% but we only allow it to be at most three (due to hard coding).
%
% The corresponding DAE Ex = Ax can have strangeness index bigger than
% zero (differentiation index bigger than one). However, note that bigger
% errors might occur if the strangeness index is too big, because
% derivatives are approximated by finite differences. We suppose that the
% strangeness index is not bigger than three.
%
% The time stepping is based on the method of steps and using Radau IIA 
% collocation.
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
%   MaxIter     Upper bound for the total number of time steps (including 
%               rejected time steps).
%   MaxReject   Upper bound for the number of rejections per time step.
%   MaxCorrect  Upper bound for the number of correction steps when using
%               long steps (step size bigger than the lag).
%
%   InitStep    Inital step size.
%   MinStep     Lower bound for the step size, default: 0.
%   MaxStep     Upper bound for the step size, default: inf.
%
%   AbsTol      Absolute tolerance, default: 1e-5.
%   RelTol      Relative tolerance, default: 1e-5.
%   LagTol      Set x(t-tau(t)):=x(t) for tau(t)<=LagTol, default: 1e-5.
%
%   StrIdx      Lower bound for the strangeness index.
%   MaxStrIdx   Upper bound for the strangeness index.
%   Shift       Lower bound for the shift index.
%   MaxShift    Upper bound for the shift index.
%
%   InitVal     Initial value, not necessarily consistent.
%
% @supporting functions:
%   timeStep
%   getRegularizedSystem
%   inflateEA
%   inflateB
%   inflatef
%   solveLinSystem
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

%-------------------------------------------------------------------------%
% set missing fields in options
%-------------------------------------------------------------------------%
if ~exist('options','var'),options = {}; end

% bounds for the itererations in the main loop
if ~isfield(options,'MaxIter')   options.MaxIter = 10000; end
if ~isfield(options,'MaxReject') options.MaxReject = 100; end
if ~isfield(options,'MaxCorrect')options.MaxCorrect = 10; end

% initial step size and bounds for the step size
if ~isfield(options,'InitStep')  options.InitStep = diff(tspan)/100; end
if ~isfield(options,'MinStep')   options.MinStep = 0; end
if ~isfield(options,'MaxStep')   options.MaxStep = inf; end

% tolerances
if ~isfield(options,'AbsTol')    options.AbsTol = 1e-5; end
if ~isfield(options,'RelTol')    options.RelTol = 1e-5; end
if ~isfield(options,'LagTol')    options.LagTol = 1e-5; end

% the DDAE's indeces (guesses, if unknown)
if ~isfield(options,'StrIdx')    options.StrIdx = 0; end
if ~isfield(options,'MaxStrIdx') options.MaxStrIdx = 3; end
if ~isfield(options,'Shift')     options.Shift = 0; end
if ~isfield(options,'MaxShift')  options.MaxShift = 3; end

% initial value (not necessarily consistent)
if ~isfield(options,'InitVal')   options.InitVal = phi(tspan(1)); end

%-------------------------------------------------------------------------%
% defining some parameters
%-------------------------------------------------------------------------%
t0 = tspan(1);
if not(isa(E,'function_handle'))
    error('E must be a function handle.'); 
end
[m,n] = size(E(0));
h = options.InitStep;
N = floor(diff(tspan)/h);
x0 = options.InitVal;

info.Rejected_Steps = 0;

%-------------------------------------------------------------------------%
% some more input checks
%-------------------------------------------------------------------------%
% checking tau
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
% checking the options
if options.MaxShift<options.Shift
    error('MaxShift must not be less than Shift.')
end
if options.MaxStrIdx<options.StrIdx
    error('MaxStrIdx must not be less than StrIdx.')
end


%-------------------------------------------------------------------------%
% preparation for the main loop
%-------------------------------------------------------------------------%
% The vector c comes from the Butcher tableau of the 3-stage Radar IIa
% method.
c=[(4-sqrt(6))/10; (4+sqrt(6))/10; 1];

% The column vector x(:,j) with length 3*n contains an approximation of the
% DDAE's solution at the time points t(j-1)+c(1)*h, t(j-1)+c(2)*h and t(j)
% in the first n, second n and third n entries, respectively, where h
% denotes the (current) step size.
x=nan(3*n,N+1);
x(1:n,1)=phi(t0+(c(1)-1)*h);
x(n+1:2*n,1)=phi(t0+(c(2)-1)*h);
t=nan(1,N+1);
t(1)=tspan(1);

%-------------------------------------------------------------------------%
% finding consistent initial value
%-------------------------------------------------------------------------%
hist_funcs = cell(l,1);
for k=1:l
    hist_funcs{k} = phi;
end
BXTAUF = @(s) B(s)*eval_all_hist_funcs(hist_funcs,s,tau(s),n,l)+f(s);
[E1,A1,~,A2,g2,mu,Z1,Z2] = regularize_strange_ldae(E,A,BXTAUF,t0,options);
x(2*n+1:3*n,1)=x0-pinv(A2)*(A2*x0+g2);
info.Strangeness_Index = mu;

%-------------------------------------------------------------------------%
% main loop - time integration
%-------------------------------------------------------------------------%
for i=1:options.MaxIter
    % Compute approximation of x at t = t(i)+h.
    
    absErr = inf;
    relErr = inf;
    
    for j=1:options.MaxReject
        % Estimate the local error by performing a full step and two half
        % steps.
        x_full = timeStep(E,A,B,f,tau,phi,t(1:i),x(:,1:i),h,options);
        x_half = timeStep(E,A,B,f,tau,phi,t(1:i),x(:,1:i),h/2,options);
        x_half = timeStep(E,A,B,f,tau,phi,[t(1:i),t(i)+h/2],[x(:,1:i),x_half],h/2,options);
        absErr = norm(x_half(2*n+1:3*n)-x_full(2*n+1:3*n));
        relErr = absErr/norm(x_half(2*n+1:3*n));
        % If the error fulfills the prescribed tolerances or the step size
        % is already equal to the minimal step size, then the step is
        % accepted. If not, the step is "rejected", the step size halved,
        % and we repeat the procedure.
        if (absErr<=options.AbsTol && relErr<=options.RelTol) || (h<=options.MinStep)
            info.Rejected_Steps = info.Rejected_Steps + j-1;
            break;
        end
        h = max(h/2,options.MinStep);
    end
    
    % Use x_half for the approximation at t(i)+c(3)*h = t(i)+h.
    x(:,i+1) = [x_full(1:2*n);x_half(2*n+1:3*n)];
    t(i+1) = t(i) + h;
    
    % Estimate the next step size h.
    h_newAbs = 0.9*h*(options.AbsTol/absErr)^(1/6);
    h_newRel = 0.9*h*(options.RelTol/relErr)^(1/6);
    h = min([h_newAbs,h_newRel,2*h]);
    
    % Impose lower and upper bounds on the step size h.
    h = max(h,options.MinStep);
    h = min([h,options.MaxStep,tspan(2)-t(i+1)]);
    
    if t(i+1)>=tspan(2)
        break
    end
end

% "cutting out" the approximate solution at t
x=x(2*n+1:3*n,1:i+1);
t=t(1:i+1);

%-------------------------------------------------------------------------%
% primary supporting functions
%-------------------------------------------------------------------------%
function x_next = timeStep(E,A,B,f,tau,phi,t,x,h,options)
% Performs EITHER a usual step of the method of steps with index reduction 
% OR a long step, i.e. h>tau and x(t-tau) has to be predicted using
% extrapolation of the last computed cubic polynomial, with index
% reduction. After extrapolating and computing the current cubic
% polynomial, we will eventually get a new x(t-tau(t)) that differs from
% the 
[m,n] = size(E(0));
% The vector c comes from the Butcher tableau of the 3-stage Radar IIa
% method.
c=[(4-sqrt(6))/10; (4+sqrt(6))/10; 1];
% Containers for the matrices of the local strangeness-free formulations at
% t_ij = T(i)+h*c(j).
Etij=nan(n,n,3);
Atij=nan(n,n,3);
Btij=nan(n,n,3);
ftij=nan(n,3);
gtij=nan(n,3);
xtau=nan(n,3);

% The FOR-loop is only used for long steps, i.e. if x(t-tau) has to be
% extrapolated. During the loop, x(t-tau) will be corrected.
for i=1:options.MaxCorrect
    if i==1
        % For each collocation point...
        isLongStep = false(); TODO
        for j=1:3
            tau_j = tau(t(end)+c(j)*h);
            if tau_j<0
                error('THE DELAY tau IN x(t-tau) IS NEGATIVE!');
            elseif tau_j<=options.LagTol;
                % If the delay is vanishing or becoming to small, then we
                % just interpret x(t-tau) as x(t).
                xtau(:,j) = zeros(n,1);
                % Calculate locally regularized form at t = t(i)-c(j)*h.
                [E1,A1,~,f1,A2,~,f2] = getRegularizedSystem(E,@(t)A(t)+B(t),@(t)zeros(m,n),f,tau,t(end)+c(j)*h,options);
                Etij(:,:,j)=[E1;zeros(size(A2))];
                Atij(:,:,j)=[A1;A2];
                Btij(:,:,j)=zeros(n,n);
                ftij(:,j)=[f1;f2];
            else
                % The delay is big enough and determine x(t-tau).
                t_tau = t(end)+c(j)*h-tau_j;
                if t_tau<t(1)
                    % t-tau(t) is less than t0, so we use the history function.
                    xtau(:,j) = phi(t_tau);
                elseif t_tau*(1+eps)<t(end)
                    % t-tau(t) is in the interval [t(1),t(end)]
                    % find the biggest time node smaller than t_i+c_j*h-tau
                    L = find(t_tau<t,1)-1;
                    % If t_i+c_j*h-tau is not a node point, i.e. not in t, then
                    % we have to interpolate. We use a polynomial of degree 3, 
                    % so we need 4 data points.
                    x0_tau=x(2*n+1:3*n,L);
                    X_tau=reshape(x(:,L+1),n,3);
                    h_tau = t(L+1)-t(L);
                    % Interpolate with Neville-Aitken.
                    xtau(:,j) = nevilleAitken(t(L)+[0;c]*h_tau,[x0_tau,X_tau],t_tau);
                else
                    % t-tau(t) is greater than t(end), use extrapolation.                    
                    %disp('Performing long step.')
                    isLongStep(j) = true;
                    if size(x,2)<2
                        warning('NOT ENOUGH POINTS FOR EXTRAPOLATION')
                        x_next=inf(3*n,1);
                        return
                    end
                    x0_tau=x(2*n+1:3*n,end-1);
                    X_tau=reshape(x(:,end),n,3);
                    h_tau = t(end)-t(end-1);
                    % Extrapolate with Neville-Aitken.
                    xtau(:,j) = nevilleAitken(t(end-1)+[0;c]*h_tau,[x0_tau,X_tau],t_tau);
                end
                % Calculate locally regularized form at t = t(i)-c(j)*h.
                [E1,A1,B1,f1,A2,B2,f2] = getRegularizedSystem(E,A,B,f,tau,t(end)+c(j)*h,options);
                Etij(:,:,j)=[E1;zeros(size(A2))];
                Atij(:,:,j)=[A1;A2];
                Btij(:,:,j)=[B1;B2];
                ftij(:,j)=[f1;f2];
            end
        end
    end
    for j=1:3
        gtij(:,j)=Btij(:,:,j)*xtau(:,j)+ftij(:,j);
    end
    % Solve the linear system.
    x_next = solveLinSystem(Etij,Atij,gtij,x,h);
    
    % No correction needed if we did not perform a long step.
    if sum(isLongStep)==0
        return
    end
    
    % Otherwise we can correct the extrapolated x(t-tau) by evaluating the
    % just computed new cubic polynomial.
    xtau_corrected = xtau;
    for j=1:3
        t_tau = t(end)+c(j)*h-tau(t(end)+c(j)*h);
        x0_tau=x(2*n+1:3*n,end);
        X_tau=reshape(x_next,n,3);
        % interpolate with Neville-Aitken
        xtau_corrected(:,j) = nevilleAitken(t(end)+[0;c]*h,[x0_tau,X_tau],t_tau);
    end
    
    % If the corrected x(t-tau) differs too much from the extrapolated one,
    % we recompute the polynomial using xtau_corrected instead of xtau. 
    if or(norm(xtau_corrected(:,isLongStep)-xtau(:,isLongStep))/norm(xtau_corrected(:,isLongStep))*h<options.RelTol,norm(xtau_corrected(:,isLongStep)-xtau(:,isLongStep))*h<options.AbsTol)
        %fprintf('Corrected after %d steps.\n',i)
        return
    else
        if i==options.MaxCorrect
        fprintf('Correction of x(t-tau(t)) failed after MaxCorrect=%d iterations.\nRemaining relative residual: %e\nRemaining absolute residual: %e\n',i,norm(xtau_corrected(:,isLongStep)-xtau(:,isLongStep))/norm(xtau_corrected(:,isLongStep))*h,norm(xtau_corrected(:,isLongStep)-xtau(:,isLongStep))*h)
        end
        xtau(:,isLongStep)=xtau_corrected(:,isLongStep);
    end
end
function [E_1,A_1,B_1,f_1,A_2,B_2,f_2,mu,K] = getRegularizedSystem(E,A,B,f,tau,ti,options)
%TODO more comments
isConst = 0;
K0 = 0;
KMax = 3;
mu0 = 0;
muMax = 3;
tolA = 1e-7;
tolR = 1e-7;
if exist('options','var')
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
            NMP((1:(mu+1)*m)+K*(muMax+1)*m,(1:(mu+1)*n)+(K-1)*(muMax+2)*n) = -inflateB(B,t_shifted(K+1),tau,mu,tolR,m,n);
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
        P = U1(1:(mu+1)*m,:)'*inflateB(B,ti,tau,mu,tolR,m,n);
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

%-------------------------------------------------------------------------%
% auxiliary functions for getRegularizedSystem()
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
function P = inflateB(B,tK,tau,mu,tol,m,n)
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

%-------------------------------------------------------------------------%
% auxiliary function for solving the linear system
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
function px = nevilleAitken(X,F,x)
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