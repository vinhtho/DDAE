function [t,x,info] = colddae(E,A,B,f,tau,phi,tspan,options)
%COLDDAE numerical solver for either
% causal causal linear delay differential-algebraic equations of the form
%   E(t)\dot{x}(t) = A(t)x(t) + \sum_{i=1}^k B_i(t)x(t-\tau_i(t)) + f(t)  for t\in(t0,tf]
%             x(t) = \phi(t),                                             for t<=t0,
% or non-causal linear delay differential-algebraic equations of the form
%   E(t)\dot{x}(t) = A(t)x(t) + B(t)x(t-\tau(t)) + f(t)  for t\in(t0,tf]
%             x(t) = \phi(t),                            for t<=t0
% with E,A,B,B_i,f,\tau,\tau_i,\phi sufficiently smooth.
%
% The strageness index can be greater than zero (see [1] for its definition
% for DDAEs). However, note that the bigger the strangeness index is, the
% larger possible numerical errors can be, especially if the all
% derivatives are approximated by finite differences. We suppose that the
% strangeness index is not greater than three due to hard coding in the 
% subroutine inflateB. If the user provides the derivative array, then the
% strangeness index can be arbitrary.
%
% The time stepping is based on the method of steps and using Radau IIA 
% collocation. The step size is chosen in each step to make the local error
% satisfy the given relative AND absolute tolerances. Also, we allow
% so-called "long steps", i.e., the step size is allowed to become bigger
% than the delay. 
%
% @parameters:
%   E,A         Coefficients of the DDAE, m-by-n matrix functions.
%   B           Coefficients of the DDAE, m-by-k*n matrix functions, where
%               k is the number of delays
%   f           m-by-1 vector function.
%   tau         Variable delay or lag,  1-by-k vector function.
%   phi         n-by-1 history function.
%   tspan       Considered time interval [t0,tf].
%   options     Struct for optional parameters, set by
%               'options.FieldName = FieldValue', see below
%
% @options
%   MaxIter     Upper bound for the total number of time steps (excluding 
%               rejected time steps).
%   MaxReject   Upper bound for the number of rejections per time step.
%   MaxCorrect  Upper bound for the number of correction steps when using
%               long steps (step size bigger than the lag), default: 10.
%
%   InitStep    Inital step size.
%   MinStep     Lower bound for the step size, default: 0.
%   MaxStep     Upper bound for the step size, default: inf.
%
%   AbsTol      Absolute tolerance, default: 1e-5.
%   RelTol      Relative tolerance, default: 1e-5.
%
%   StrIdx      Lower bound for the strangeness index, default: 0.
%   MaxStrIdx   Upper bound for the strangeness index, default: 3
%   Shift       Lower bound for the shift index, default: 0.
%   MaxShift    Upper bound for the shift index, default: 3.
%
%   InitVal     Initial value, not necessarily consistent, default: 
%               phi(tspan(1)).
%
%   IsConst     A boolean, true if E,A,B,tau are constant (then the
%               regularized system is computed only once, i.e. the
%               solver needs less computation time), default: false.
%   
%   DArray      A struct containing M, P, and g as function handles in the
%               first, second and third entry, resp, which enables the user
%               to pass the exact derivative array of the DDAE, which is of
%               the form M(t)*z(t)=P(t)*z(t-tau(t))+g(t), whith
%               z=[x;\dot{x};\ddot{x};...;x^{(mu+1)}], default: not set.
%               
% @supporting functions:
%   timeStep
%   solveLinSystem
%   nevilleAitken
%   newton
%   getRegularizedSystem
%   inflateEA
%   inflateB
%   inflateAndShiftf
%   matrixDifferential
%   null2
%   orth2
%
% @return values:
%   t           t(i+1) = t0+h_i with h_i the i-th step size.
%   x           numerical solution at the time nodes in t.
%   info        Struct with information.
%
% @author:
%       Vinh Tho Ma, TU Berlin, mavinh@math.tu-berlin.de
%       Phi Ha, TU Berlin, ha@math.tu-berlin.de
%
% References:
%  [1]  Phi Ha, Vinh Tho Ma and Volker Mehrmann: A new solver for linear
%       delay differential algebraic equations, Technical Report, Institut
%       fuer Mathematik, TU Berlin, Berlin, Germany, 2015.

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

% the DDAE's indeces (guesses, if unknown)
if ~isfield(options,'StrIdx')    options.StrIdx = 0; end
if ~isfield(options,'MaxStrIdx') options.MaxStrIdx = 3; end
if ~isfield(options,'Shift')     options.Shift = 0; end
if ~isfield(options,'MaxShift')  options.MaxShift = 3; end

% initial value (not necessarily consistent)
if ~isfield(options,'InitVal')   options.InitVal = phi(tspan(1)); end

% Are E, A, B, and tau constant functions?
if ~isfield(options,'IsConst')   options.IsConst = false; end

% Checking the d-array
options.DArray_provided = 0;
if isfield(options,'DArray')
    if numel(options.DArray) == 3
        options.DArray_provided = 1;
    end
end

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

% predefining info's fields
info.Strangeness_index = -1;
info.Shift_index = -1;
info.Number_of_differential_eqs = -1;
info.Number_of_algebraic_eqs = -1;
info.Rejected_steps = 0;
info.Computation_time = -1;

%-------------------------------------------------------------------------%
% some more input checks
%-------------------------------------------------------------------------%
% checking tau
if not(isa(tau,'function_handle'))
    error('Delay tau must be a function handle.');
end
k = numel(tau(tspan(1)));
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
if or(size(B(0),1)~=m,size(B(0),2)~=k*n)
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
[E1,A1,B1,~,A2,B2,f2,mu,K,U,Z1,Z2] = getRegularizedSystem(E,A,B,f,tau,t0,options,m,n);
xtau0 = nan(k*n,1);
tau0 = tau(0);
for l=1:k
    xtau0((l-1)*n+1:l*n)=phi(t0-tau0(l));
end
x(2*n+1:3*n,1)=x0-pinv(A2)*(A2*x0+B2*xtau0+f2);

% Store the just computed regularized system for later use if E,A,B,tau are
% constant.
if options.IsConst
    options.RegSystem = {E1,A1,B1,A2,B2,mu,K,U,Z1,Z2};
end

% Initialize info's fields and set them by using the just computed regular
% system.
info.Strangeness_index = mu;
info.Shift_index = K;
info.Number_of_differential_eqs = size(E1,1);
info.Number_of_algebraic_eqs = size(A2,1);

%-------------------------------------------------------------------------%
% Main loop.
%-------------------------------------------------------------------------%
tic
for i=1:options.MaxIter
    % Compute approximation of x at t = t(i)+h.
    
    % The absolute and relative tolerances on the local error.
    absErr = inf;
    relErr = inf;
    
    for j=1:options.MaxReject
        % Estimate the local error by performing a full step and two half
        % steps.
        [x_full,info] = timeStep(E,A,B,f,tau,phi,t(1:i),x(:,1:i),h,options,m,n,info);
        x_half = timeStep(E,A,B,f,tau,phi,t(1:i),x(:,1:i),h/2,options,m,n,info);
        x_half = timeStep(E,A,B,f,tau,phi,[t(1:i),t(i)+h/2],[x(:,1:i),x_half],h/2,options,m,n,info);
        absErr = norm(x_half(2*n+1:3*n)-x_full(2*n+1:3*n));
        relErr = absErr/norm(x_half(2*n+1:3*n));
        
        % If the error fulfills the prescribed tolerances or the step size
        % is already equal to the minimal step size, then the step is
        % accepted. If not, the step is "rejected", the step size halved,
        % and we repeat the procedure.
        if (absErr<=options.AbsTol && relErr<=options.RelTol) || (h<=options.MinStep)
            info.Rejected_steps = info.Rejected_steps + j-1;
            break;
        end
        h = max(h/2,options.MinStep);
    end
    
    % Use x_half for the approximation at t(i)+c(3)*h = t(i)+h and the
    % other 2*n entries from x_full (could also set x(:,i+1)=x_full in
    % principle).
    x(:,i+1) = [x_full(1:2*n);x_half(2*n+1:3*n)];
    t(i+1) = t(i) + h;
    
    % Estimate the next step size h, but don't let it grow more than 2*h.
    %
    % TODO: Wrong formula at the moment. The error of the delay terms are
    % not considered yet.
    %
    h_newAbs = 0.9*h*(options.AbsTol/absErr)^(1/6);
    h_newRel = 0.9*h*(options.RelTol/relErr)^(1/6);
    h = min([h_newAbs,h_newRel,2*h]);
    
    % Impose lower and upper bounds on the step size h, if any.
    h = max(h,options.MinStep);
    h = min([h,options.MaxStep,tspan(2)-t(i+1)]);
    
    if t(i+1)>=tspan(2)
        break
    end
end

% "Cutting out" the approximate solution at the time points given in t,
% i.e. t0+c(3)*h_i, i = 1,2,....
x=x(2*n+1:3*n,1:i+1);
t=t(1:i+1);
info.Computation_time = toc;

%-------------------------------------------------------------------------%
% Supporting functions.
%-------------------------------------------------------------------------%
function [x_next,info] = timeStep(E,A,B,f,tau,phi,t,x,h,options,m,n,info)
% Performs EITHER a usual step of the method of steps with index reduction 
% OR a long step, i.e. h>tau and x(t-tau) has to be predicted using
% extrapolation of the last computed cubic polynomial, with index
% reduction. After extrapolating and computing the current cubic
% polynomial, we will eventually get a new x(t-tau(t)) by interpolating 
% that differs from the extrapolated x(t-tau(t)). We use the new computed
% x(t-tau(t)) to compute a hopefully more exact cubic polynomial. This can
% be done over and over again and hopefully this process converges.

%-------------------------------------------------------------------------%
% Preparations
%-------------------------------------------------------------------------%
% The vector c comes from the Butcher tableau of the 3-stage Radar IIa
% method.
c=[(4-sqrt(6))/10; (4+sqrt(6))/10; 1];
% Containers for the matrices of the local strangeness-free formulations at
% t_ij = T(i)+h*c(j).
Etij=nan(n,n,3);
Atij=nan(n,n,3);
k = numel(tau(t(1)));
Btij=nan(n,k*n,3);
ftij=nan(n,3);
btij=nan(n,3);
xtau=nan(k*n,3);
DArray_provided = options.DArray_provided;


% A matrix of booleans used for logical indexing, used in order to
% adress only those entries which have been computed by long steps.
isLongStep = false(k*n,3);

%-------------------------------------------------------------------------%
% Main loop.
%-------------------------------------------------------------------------%
% The FOR-loop is only used for long steps, i.e. if x(t-tau) has to be
% extrapolated. During the loop, x(t-tau) will be corrected.
for i=1:options.MaxCorrect
    if i==1
        % j iterates through the stages, i.e. t(end) + c(j)*h.
        for j=1:3
            % Calculate locally regularized form at t=t(end)-c(j)*h. If
            % E,A,B,tau are constant, then we use the regularized system
            % that we computed before during the computation of the
            % consistent initial value. Otherwise compute the locally 
            % regularized system again by calling getRegularizedSystem.
            if options.IsConst
                E1 = options.RegSystem{1};
                A1 = options.RegSystem{2};
                B1 = options.RegSystem{3};
                A2 = options.RegSystem{4};
                B2 = options.RegSystem{5};
                mu = options.RegSystem{6};
                K = options.RegSystem{7};
                U = options.RegSystem{8};
                Z1 = options.RegSystem{9};
                Z2 = options.RegSystem{10};
                tau_const = tau(t(1));
                if DArray_provided
                    % Use provided derivative array, if provided and enough
                    % derivatives.
                    g = nan((mu+1)*(K+1)*m,1);
                    for p = 0:K
                        g_provided = feval(options.DArray{3},t(end)+c(j)*h+p*tau_const);
                        if length(g_provided)<mu
                            warning('PLEASE PROVIDE MORE DERIVATIVES OR USE NO DERIVATIVE ARRAYS.')
                            g = inflateAndShiftf(f,t(end)+c(j)*h+(0:K)*tau_const,K,mu,options.RelTol,m);
                            break;
                        end
                        g((1:(mu+1)*m)+p*(mu+1)*m) = g_provided(1:(mu+1)*m);
                    end
                else
                    % Else use finite differences.
                    g = inflateAndShiftf(f,t(end)+c(j)*h+(0:K)*tau_const,K,mu,options.RelTol,m);
                end
                f1 = Z1'*U'*g;
                f2 = Z2'*U'*g;
            else
                [E1,A1,B1,f1,A2,B2,f2,mu,K] = getRegularizedSystem(E,A,B,f,tau,t(end)+c(j)*h,options,m,n);
                if mu<info.Strangeness_index,info.Strangeness_index=mu;end
                if K<info.Shift_index,info.Strangeness_index=K;end
            end
            
            % Store the locally regularized system at t = t(end)+c(j)+h.
            Etij(:,:,j)=[E1;zeros(size(A2))];
            Atij(:,:,j)=[A1;A2];
            Btij(:,:,j)=[B1;B2];
            ftij(:,j)=[f1;f2];
            
            % Now compute x(t-tau_1),...,x(t-tau_k) either by interpolating
            % or by extrapolating.
            tau_j = tau(t(end)+c(j)*h);
            % k iterates through all delays tau_k(t(end)+c(j)*h).
            for l = 1:k
                if tau_j(l)<=0
                    error('THE DELAY IS NOT POSITIVE!');
                else
                    t_tau = t(end)+c(j)*h-tau_j(l);
                    if t_tau<t(1)
                        % t-tau(t) is less than t0, so we use the history
                        % function.
                        xtau((l-1)*n+1:l*n,j) = phi(t_tau);
                    elseif t_tau*(1+eps)<t(end)
                        % t-tau(t) is in the interval [t(1),t(end)]. Find
                        % the biggest time node smaller than t_i+c_j*h-tau
                        L = find(t_tau<t,1)-1;
                        
                        % If t(end)i+c(j)*h-tau is not a node point, i.e.
                        % not in t, then we have to interpolate. We use a
                        % polynomial of degree 3, so we need 4 data pairs.
                        x0_tau = x(2*n+1:3*n,L);
                        X_tau = reshape(x(:,L+1),n,3);
                        h_tau = t(L+1)-t(L);
                        % Interpolate with Neville-Aitken's algorithm.
                        xtau((l-1)*n+1:l*n,j) = nevilleAitken(t(L)+[0;c]*h_tau,[x0_tau,X_tau],t_tau);
                    else
                        % t-tau(t) is greater than t(end), so we use 
                        % extrapolation.
                        isLongStep((l-1)*n+1:l*n,j) = true;
                        if size(x,2)<2
                            % We don't have enough extrapolation points,
                            % this happens if the solver wants to perform a
                            % long step right at the beginning, when no
                            % cubic polynomial has been computed yet.
                            x_next=inf(3*n,1);
                            return
                        end
                    x0_tau=x(2*n+1:3*n,end-1);
                    X_tau=reshape(x(:,end),n,3);
                    h_tau = t(end)-t(end-1);
                    % Extrapolate with Neville-Aitken's algorithm.
                    xtau((l-1)*n+1:l*n,j) = nevilleAitken(t(end-1)+[0;c]*h_tau,[x0_tau,X_tau],t_tau);
                    end
                end
            end
        end
    end
    
    % Compute b(t) := B(t)*x(t-tau)+f(t) at t = t(end)+c(j)*h, j = 1,2,3.
    for j=1:3
        btij(:,j)=Btij(:,:,j)*xtau(:,j)+ftij(:,j);
    end
    
    % Solve the linear system.
    x_next = solveLinSystem(Etij,Atij,btij,x(:,end),h);
    
    % No correction needed if we did not perform a long step.
    if sum(sum(isLongStep))==0
        return
    end
    
    % Otherwise we can correct the extrapolated x(t-tau) by evaluating the
    % just computed cubic polynomial given in x_next.
    xtau_corrected = xtau;
    for j=1:3
        % Note that t_tau is a vector in the multiple delay case!
        t_tau = t(end)+c(j)*h-tau(t(end)+c(j)*h);
        x0_tau=x(2*n+1:3*n,end);
        X_tau=reshape(x_next,n,3);
        % Interpolate with Neville-Aitken, if x(t-tau_k) was computed with
        % a long step.
        for l=1:k
            if isLongStep(l*n,j)
                xtau_corrected((l-1)*n+1:l*n,j) = nevilleAitken(t(end)+[0;c]*h,[x0_tau,X_tau],t_tau(l));
            end
        end
    end
    
    % If xtau_corrected differs too much from the xtau, then
    % we recompute the polynomial using xtau_corrected. 
    if or(norm(xtau_corrected(isLongStep)-xtau(isLongStep))/norm(xtau_corrected(isLongStep))*h<options.RelTol,norm(xtau_corrected(isLongStep)-xtau(isLongStep))*h<options.AbsTol)
        %fprintf('Corrected after %d steps.\n',i)
        return
    else
        if i==options.MaxCorrect
        fprintf('Correction of x(t-tau(t)) failed after MaxCorrect=%d iterations.\nRemaining relative residual: %e\nRemaining absolute residual: %e\n',i,norm(xtau_corrected(isLongStep)-xtau(isLongStep))/norm(xtau_corrected(isLongStep))*h,norm(xtau_corrected(isLongStep)-xtau(isLongStep))*h)
        end
        xtau(isLongStep)=xtau_corrected(isLongStep);
    end
end

function x_next = solveLinSystem(Etij,Atij,ftij,xi,h)
% Supporting function for solving the linear system to determine the next
% cubic polynomial, which approximates the solution at some time points of
% the form t+c(j)*h, j=1,2,3.
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
    bb((j-1)*n+1:j*n)=ftij(:,j)-Etij(:,:,j)/h*v0(j)*xi(2*n+1:3*n);
end
% the solution is a vector with length 3*n, it consists of the 3 values
% of the polynomial at the collocation points t(i)+c(j)*h, j=1..3
x_next=AA\bb;

function px = nevilleAitken(X,F,x)
% Evaluate the (vector valued) polynomial given by (X,F) at the time point 
% x.
n=length(X);
px=F;
for i=1:n-1
    for j=1:n-i
        px(:,j)=((x-X(j))*px(:,j+1)-(x-X(j+i))*px(:,j))/(X(j+i)-X(j));
    end
end
px=px(:,1);

function [E_1,A_1,B_1,f_1,A_2,B_2,f_2,mu,K,U,Z1,Z2] = getRegularizedSystem(E,A,B,f,tau,ti,options,m,n)
% Calculates the regular system of the DDAE given by the triple (E,A,B)
% locally at the point t=ti.

%-------------------------------------------------------------------------%
% Preparation.
%-------------------------------------------------------------------------%
KMax=options.MaxShift;
muMax=options.MaxStrIdx;
% We use AbsTol for calculating null spaces, ranges and use it in the
% Newton method.
AbsTol=options.AbsTol; 
% RelTol is used for matrix_differential
RelTol=options.RelTol; 
K0=max(0,options.Shift);
mu0=options.StrIdx;
k = numel(tau(ti));
muMax = max(mu0,muMax);
DArray_provided = options.DArray_provided;


% Container for the Double Inflated SYStem. If the user does not know the
% exact strangeness index or the shift index, then the solver checks all
% combinations from, say, mu = 0,1,2,3 and K = 0,1,2,3, which corresponds
% to the iterative computation of 9 different double inflated systems.
% However, computing every double inflated system from scratch would lead
% to recomputations of already computed values, e.g. the system for K=2 and
% mu = 0 contains 2*m equations of the system for K=1 and mu=3. That is why
% everything is stored in DISYS, and appropriate entries are selected by
% the logical indeces idxN, idxM, idx1 and idx2.
DISYS = zeros((muMax+1)*m*(KMax+1),(muMax+2)*n*(KMax+1));

% Logical indices of to the coefficients of x and \dot{x} in DISYS.
idx2M = false((muMax+2)*n*(KMax+1),1);
idx2M(n+1:2*n) = true;
idx2N = false((muMax+2)*n*(KMax+1),1);
idx2N(1:n) = true;
% Compute t(K+1) = t' which fulfills t'-tau(t') = t(K). Works
% only for single delay.
t_shifted = nan(KMax+1,1);
t_shifted(1) = ti;
if k==1
%     for K=1:KMax
%         t_shifted(K+1) = newton(@(t1) t1-tau(t1)-t_shifted(K),t_shifted(K)+tau(t_shifted(K)),1e-15,10);
%     end
    for i = 1:KMax
        t_shifted(i+1) = fsolve(@(t1) t1-tau(t1)-t_shifted(i),t_shifted(i)+tau(t_shifted(i)),optimset('Display','off'));
    end    
end

%-------------------------------------------------------------------------%
% Main loops: Increasing the shift index K and strangeness index mu until
% we get a regular system.
%-------------------------------------------------------------------------%
for K = 0:KMax
    %---------------------------------------------------------------------%
    % Compute the (current) double inflated system.
    %---------------------------------------------------------------------%
    % Logical indices for selecting certain rows and columns in DISYS.
    idx1 = false((muMax+1)*m*(KMax+1),1);
    idx2 = false((muMax+2)*n*(KMax+1),1);
    
    mu_provided = inf;
    % Use user provided derivative array, if provided.
    if DArray_provided
        inflateEA_provided = feval(options.DArray{1},t_shifted(K+1));
        inflateB_provided = feval(options.DArray{2},t_shifted(K+1));
        mu_provided = size(inflateEA_provided,1)/m-1;
    end
    
    for mu = 0:muMax
        if mu>mu_provided
            warning('PLEASE PROVIDE MORE DERIVATIVES OR USE NO DERIVATIVE ARRAYS.');
        end
        if DArray_provided && mu<=mu_provided
            DISYS((1:(mu+1)*m)+K*(muMax+1)*m,(1:(mu+2)*n)+K*(muMax+2)*n) = inflateEA_provided((1:(mu+1)*m),(1:(mu+2)*n));
        else
            DISYS((1:(mu+1)*m)+K*(muMax+1)*m,(1:(mu+2)*n)+K*(muMax+2)*n) = inflateEA(E,A,t_shifted(K+1),mu,RelTol);
        end
        if K>0
            if k>1
                error('REGULARIZATION FOR NONCAUSAL DDAES WITH MULTIPLE DELAYS NOT IMPLEMENTED YET.')
            end
            if isfield(options,'DArray') && mu<=mu_provided
                DISYS((1:(mu+1)*m)+K*(muMax+1)*m,(1:(mu+1)*n)+(K-1)*(muMax+2)*n) = -inflateB_provided(1:(mu+1)*m,1:(mu+1)*n);
            else
                DISYS((1:(mu+1)*m)+K*(muMax+1)*m,(1:(mu+1)*n)+(K-1)*(muMax+2)*n) = -inflateB(B,t_shifted(K+1),tau,mu,RelTol,m,n);
            end
            if K<K0
                continue;
            end
            % logical indices for the coefficients of x(t+tau), x(t+2*tau),
            % ..., x(t+K*tau) and its derivatives up to order mu in DISYS
            for i = 1:K
                idx1((1:(mu+1)*m)+i*(muMax+1)*m)=true;
                idx2((1:(mu+2)*n)+i*(muMax+2)*n)=true;
            end
        end
        if mu<mu0
            continue;
        end
        % logical indices for the coefficients of x and its derivatives up 
        % to order mu in DISYS
        idx1(1:(mu+1)*m) = true;
        idx2((2*n+1):(mu+2)*n) = true;
        
        % Extract a system without x(t+tau), x(t+2*tau), ..., x(t+K*tau)
        % and its derivatives.
        U1 = null2(DISYS(idx1,idx2)',AbsTol);
        % Apperently, the two lines below lead to an error, the intention
        % was to eliminate redundancy, but P is not considered here, so
%         U2 = orth2(U1'*DISYS(idx1,or(idx2N,idx2M)),tol);
%         U=U1*U2;
        M = U1'*DISYS(idx1,idx2M);
        N = -U1'*DISYS(idx1,idx2N);
        
        % We determine the characteristic values a and d and check wether
        % we have a regular system, i.e. a+d=n.
        Z2 = null2(M',AbsTol);
        A_2 = Z2'*N;
        T2 = null2(A_2,AbsTol);
        Z1 = orth2(M*T2,AbsTol);
        E_1 = Z1'*M;
        a = rank(A_2,AbsTol);
        d = rank(E_1,AbsTol);
        if a+d<n
            continue
        end
        
        %-----------------------------------------------------------------%
        % Computation a regular system, which contains no derivatives of
        % x(t-tau) and no redundant equations.
        %-----------------------------------------------------------------%
        % We eliminate all equations which contain derivatives of
        % x(t-tau).
        if DArray_provided && mu<=mu_provided
            P = U1(1:(mu+1)*m,:)'*inflateB_provided(1:(mu+1)*m,1:(mu+1)*k*n);
        else
            P = U1(1:(mu+1)*m,:)'*inflateB(B,ti,tau,mu,RelTol,m,n);
        end
        U2 = null2(P(:,k*n+1:end)',AbsTol);
        P = U2'*P(:,1:k*n);
        M = U2'*M;
        N = U2'*N;
        
        % If the system is not advanced, then we still should be able to
        % extract a regular system with a+d=n, if not, then we very likely
        % have an advanced system.
        Z2 = null2(M',AbsTol);
        A_2 = Z2'*N;
        T2 = null2(A_2,AbsTol);
        Z1 = orth2(M*T2,AbsTol);
        E_1 = Z1'*M;
        if or(rank(A_2,AbsTol)~=a,rank(E_1,AbsTol)~=d)
            error('ACCORDING TO THE CHOSEN TOLERANCE, THE DDAE IS ADVANCED. THIS SOLVER CAN NOT HANDLE ADVANCED DDAES YET.')
        end
        
        % Remove redundant algebraic and differential equations.
        if size(A_2,1)>a
            Y2 = orth2(A_2,AbsTol);
            A_2 = Y2'*A_2;
            % update the selector Z2
            Z2 = Z2*Y2;
        end
        if size(E_1,1)>d
            Y1 = orth2(E_1,tol);
            E_1 = Y1'*E_1;
            % update the selector Z1
            Z2 = Z2*Y2;
        end
        
        % Compute the rest of the regular system.
        g = nan((mu+1)*m*(K+1),1);
        if DArray_provided && mu<=mu_provided
            for j = 0:K
                g(j*(mu+1)*m+1:(j+1)*(mu+1)*m) = feval(options.DArray{3},t_shifted(j+1));
            end
        else
            g = inflateAndShiftf(f,t_shifted,K,mu,RelTol,m);     
        end
        U=U1*U2;
        g = U'*g;
        B_2 = Z2'*P;
        f_2 = Z2'*g;
        A_1 = Z1'*N;
        B_1 = Z1'*P;
        f_1 = Z1'*g;
        return
    end
end
error('MAXIMAL NUMBER OF SHIFTS AND STRANGENESS REDUCTION STEPS REACHED. REGULARIZATION OF THE DDAE FAILED.')

function xn = newton(f,x0,AbsTol,MaxIter)
% Newton's method for finding the root of a scalar function f.
for i=1:MaxIter
    f0 = f(x0);
    if f0<=AbsTol
        xn = x0;
        return
    end
    Df0 = matrixDifferential(f,x0,1);
    xn = x0-f0/Df0;
    if abs(xn-x0)<=AbsTol
        return
    end
    x0 = xn;
end

function NM = inflateEA( E,A,t,mu,tol )
% Computes the derivative array of (E,A) by differentiating mu times.
%   INPUT
%   -----
%   E       fcn_handle      m-by-n leading matrix function
%   A       fcn_handle      m-by-n matrix function
%   t       double          the time
%   mu      double          the strangeness index
%   tol    double          the relative tolerance
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
    % Make dE and dA contain all derivatives up to order l.
    dE(l*m+1:(l+1)*m,1:n) = matrixDifferential( E,t,l);
    dA(l*m+1:(l+1)*m,1:n) = matrixDifferential( A,t,l);
    %Expand M_(l-1) to M_l.
    for j = 0:l-1
        k = l-j;
        NM(l*m+1:(l+1)*m,(j+1)*n+1:(j+2)*n) = nchoosek(l,j)*dE(k*m+1:(k+1)*m,:)-nchoosek(l,j+1) * dA((k-1)*m+1:k*m,:);
    end
    NM(l*m+1:(l+1)*m,(l+1)*n+1:(l+2)*n) = dE(1:m,1:n);
end
NM(:,1:n) = -dA;

function P  = inflateB(B,tK,tau,mu,tol,m,n)
% Builds the matrix
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
% i.e. the derivative array of B(t)*x(t-tau(t)).
%
% it is hard coded up to mu = 3 because we haven't found a nice forumla yet
% mu could be arbitrary in principle
k=numel(tau(tK));
B0 = B(tK);
switch mu
    case 0
        P = B0;
    case 1
        B1 = matrixDifferential(B,tK,1);
        tau1 = matrixDifferential(tau,tK,1);
        
        P = [
            B0,zeros(m,k*n);
            B1,B0.*kron((1-tau1),ones(m,n))
            ];
    case 2
        B1 = matrixDifferential(B,tK,1);
        B2 = matrixDifferential(B,tK,2);
        tau1 = matrixDifferential(tau,tK,1);
        tau2 = matrixDifferential(tau,tK,2);
        P = [
            B0,zeros(m,2*k*n);
            B1,B0.*kron(1-tau1,ones(m,n)),zeros(m,k*n);
            B2,2*B1.*kron(1-tau1,ones(m,n))-B0.*kron(tau2,ones(m,n)),B0.*kron((1-tau1).^2,ones(m,n))
            ];
    case 3
        B1 = matrixDifferential(B,tK,1);
        B2 = matrixDifferential(B,tK,2);
        B3 = matrixDifferential(B,tK,3);
        tau1 = matrixDifferential(tau,tK,1);
        tau2 = matrixDifferential(tau,tK,2);
        tau3 = matrixDifferential(tau,tK,3);
        P = [
            B0,zeros(m,3*k*n);
            B1,B0.*kron(1-tau1,ones(m,n)),zeros(m,2*k*n);
            B2,2*B1.*kron(1-tau1,ones(m,n))-B0.*kron(tau2,ones(m,n)),B0.*kron((1-tau1).^2,ones(m,n)),zeros(m,k*n);
            B3,3*B2.*kron(1-tau1,ones(m,n))-3*B1.*kron(tau2,ones(m,n))-B0.*kron(tau3,ones(m,n)),3*B1.*kron((1-tau1).^2,ones(m,n))-3*B0.*kron((1-tau1).*tau2,ones(m,n)),B0.*kron((1-tau1).^3,ones(m,n))
            ];
end

function g  = inflateAndShiftf(f,t_shifted,K,mu,tol,m)
% Builds the vector
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
%   |_                  _|.

g = zeros((K+1)*(mu+1)*m,1);
for j = 0:K
    for i = 1:(mu+1)
        g(((i-1)*m+1:i*m)+j*(mu+1)*m) = matrixDifferential(f,t_shifted(j+1),i-1);
    end
end

function dA = matrixDifferential(A,t,k)
% Approximates the k-th time derivative of the (matrix) function A at the
% time point t.
eps=0.01;
delta=sqrt(eps)*max(0.01,abs(t)); % TODO find a better delta which avoids cancellation
summe = A(t+(k/2)*delta);
for i=1:k
    summe = summe + (-1)^i*round(prod(((k-i+1):k)./(1:i)))*A(t+(k/2-i)*delta);
end
dA=summe/delta^k;

function Z = null2(A,AbsTol)
% Slight modification of MATLAB's null function. Singular values smaller
% than abstol are considered to be zero.
[m,n] = size(A);
[~,S,V]=svd(A,0);
if m > 1
    s = diag(S);
elseif m == 1
    s = S(1);
else s = 0;
end
r = sum(s > AbsTol);
Z = V(:,r+1:n);

function Q = orth2(A,AbsTol)
% Slight modification of MATLAB's orth function. Singular values smaller
% than abstol are considered to be zero.
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
r = sum(s > AbsTol);
Q = U(:,1:r);