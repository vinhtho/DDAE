function varargout = colddae(E,A,B,f,tau,phi,tspan,options)
% COLDDAE a collocation based solver for linear DDAEs
%
% This solver comprises of two solvers


if ~exist('options','var')
    options.IsCausal = input('Is the DDAE causal?\n 1: Yes\n 0: No\n'); 
elseif ~isfield(options,'IsCausal')
    options.IsCausal = input('Is the DDAE causal?\n 1: Yes\n 0: No\n');
end

if options.IsCausal==1
    [t,x,info]=solve_causal_ddae(E,A,B,f,tau,phi,tspan,options);
else
    [t,x,info]=solve_noncausal_ddae(E,A,B,f,tau,phi,tspan,options);
end

if nargout<=1
    sol.t = t;
    sol.x = x;
    sol.info = info;
    if nargout==0
        plot(t,x)
    end
    varargout = {sol};
else
    varargout={t,x,info};
end