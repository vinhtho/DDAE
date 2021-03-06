function [t,x] = solve_varshifted_lddae(E,A,B,f,tau,phi,tspan,options)
%solve_varshifted_lddae
%
%   Solver for linear delay differential algebraic equations (DDAEs) with
%   variable coefficients and single, variable delay tau(t) >= 0, i.e.
%           .
%       E(t)x(t) = A(t)x(t) + B(t)x(t-tau(t)) + f(t),   t0 <= t <= tf,
%           x(t) = phi(t),                              t <= t0,
%
%   based on using the method of steps and using Radau IIA collocation.
%                          .
%   The corresponding DAE Ex = Ax can have strangeness index bigger than
%   zero (differentiation index bigger than one). However, note that bigger
%   errors might occur if the strangeness index is too big, because
%   derivatives are approximated by finite differences. We suppose that the
%   strangeness index is not bigger than three.
%
%   Furthermore, we allow the so called shift index to be bigger than zero,
%   but not bigger than three (due to hard coding).
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
% set missing fields in options
%-------------------------------------------------------------------------%
if ~exist('options','var'),options = {}; end
if ~isfield(options,'AbsTol')    options.AbsTol = 1e-7; end
if ~isfield(options,'MaxIter')   options.MaxIter = 10000; end
if ~isfield(options,'MaxShift')  options.MaxShift = 3; end
if ~isfield(options,'MaxStrIdx') options.MaxStrIdx = 3; end
if ~isfield(options,'RelTol')    options.RelTol = 1e-7; end
if ~isfield(options,'Shift')     options.Shift = 0; end
if ~isfield(options,'StepSize')  options.StepSize = diff(tspan)/100; end
if ~isfield(options,'StrIdx')    options.StrIdx = 0; end
if ~isfield(options,'x0')        options.x0 = phi(tspan(1)); end

%-------------------------------------------------------------------------%
% defining some parameters
%-------------------------------------------------------------------------%
t0 = tspan(1);
E0 = E(tspan(1));
n = size(E0,2);
h = options.StepSize;
N = floor(diff(tspan)/h);
x0 = options.x0;

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
[~,~,~,~,A2,B2,f2,~,~,~,~,~] = regularize_varshift_strange_lddae(E,A,B,f,tau,t0,options);
x(2*n+1:3*n,1)=x0-pinv(A2)*(A2*x0+B2*phi(t0-tau(t0))+f2);

%-------------------------------------------------------------------------%
% main loop - time integration
%-------------------------------------------------------------------------%

for i=1:N
    % Compute approximation of x at t = t(i)+h.
    
    
    x_full = performLongStep(E,A,B,f,tau,phi,t(1:i),x(:,1:i),h,options);
    x_half = performLongStep(E,A,B,f,tau,phi,t(1:i),x(:,1:i),h/2,options);
    x_half = performLongStep(E,A,B,f,tau,phi,[t(1:i),t(i)+h/2],[x(:,1:i),x_half],h/2,options);
    
    x(:,i+1) = x_full;
    t(i+1) = t(i) + h;
    
    err = norm(x_half(2*n+1:3*n)-x_full(2*n+1:3*n));
    h = max(min([(options.RelTol/err)^(1/6),2*h,tspan(2)-t(i+1)]),h/2); 
    
    if t(i+1)>=tspan(2)
        break
    end
end

% "cutting out" the approximate solution at t
x=x(2*n+1:3*n,1:i+1);
t=t(1:i+1);
%-------------------------------------------------------------------------%
% END OF solve_varshifted_lddae
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% perform one step of the method of steps, h is smaller than tau
%
% NOTE: Might get replaced by performLongStep
%-------------------------------------------------------------------------%
function x_next = performStep(E,A,B,f,tau,phi,t,x,h,options)
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
% perform EITHER a usual step of the method of steps with index reduction 
% OR a long step, i.e. h>tau and x(t-tau) has to be predicted using
% extrapolation of the last computed cubic polynomial, with index reduction
%
% NOTE: might be renamed to performStep later
%-------------------------------------------------------------------------%
function x_next = performLongStep(E,A,B,f,tau,phi,t,x,h,options)
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
for i=1:10
    if i==1
        % For each collocation point...
        isLongStep = [false,false,false];
        for j=1:3
            tau_j = tau(t(end)+c(j)*h);
            if tau_j<0
                error('THE DELAY tau IN x(t-tau) IS NEGATIVE!');
            elseif abs(tau_j)<options.AbsTol
                % If the delay is vanishing or becoming to small, then we
                % just interpret x(t-tau) as x(t).
                xtau(:,j) = zeros(n,1);
                % Calculate locally regularized form at t = t(i)-c(j)*h.
                [E1,A1,~,f1,A2,~,f2] = regularize_varshift_strange_lddae(E,@(t)A(t)+B(t),@(t)zeros(m,n),f,tau,t(end)+c(j)*h,options);
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
                    xtau(:,j) = poleval_neville_aitken(t(L)+[0;c]*h_tau,[x0_tau,X_tau],t_tau);
                else
                    % t-tau(t) is greater than t(end), use extrapolation.                    
                    %disp('Performing long step.')
                    isLongStep(j) = true;
                    if size(x,2)<2
                        error('NOT ENOUGH POINTS FOR EXTRAPOLATION')
                    end
                    x0_tau=x(2*n+1:3*n,end-1);
                    X_tau=reshape(x(:,end),n,3);
                    h_tau = t(end)-t(end-1);
                    % Extrapolate with Neville-Aitken.
                    xtau(:,j) = poleval_neville_aitken(t(end-1)+[0;c]*h_tau,[x0_tau,X_tau],t_tau);
                end
                % Calculate locally regularized form at t = t(i)-c(j)*h.
                [E1,A1,B1,f1,A2,B2,f2] = regularize_varshift_strange_lddae(E,A,B,f,tau,t(end)+c(j)*h,options);
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
        xtau_corrected(:,j) = poleval_neville_aitken(t(end)+[0;c]*h,[x0_tau,X_tau],t_tau);
    end
    
    % If the corrected x(t-tau) differs too much from the extrapolated one,
    % we recompute the polynomial using xtau_corrected instead of xtau. 
    if or(norm(xtau_corrected(:,isLongStep)-xtau(:,isLongStep))/norm(xtau_corrected(:,isLongStep))*h<options.RelTol,norm(xtau_corrected(:,isLongStep)-xtau(:,isLongStep))*h<options.AbsTol)
        %fprintf('Corrected after %d steps.\n',i)
        return
    else
        if i==10
        fprintf('%d-th correction of x(t-tau(t)).\nRemaining relative residual: %e\nRemaining absolute residual: %e\n',i,norm(xtau_corrected(:,isLongStep)-xtau(:,isLongStep))/norm(xtau_corrected(:,isLongStep))*h,norm(xtau_corrected(:,isLongStep)-xtau(:,isLongStep))*h)
        end
        xtau(:,isLongStep)=xtau_corrected(:,isLongStep);
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