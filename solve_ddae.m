function [t,x] = solve_ddae(data,tspan,options)
if ~exist('options','var')
    options = {};
end
disp('Please choose the type of the DDAE.')
choice = input('(1) Shift index zero\n(2) Shift index bigger than zero\n');
switch choice
    case 1
        choice = input('(1) Linear\n(2) Nonlinear\n');
        switch choice
            case 1
                E = data{1};
                A = data{2};
                B = data{3};
                f = data{4};
                tau = data{5};
                phi = data{6};
                choice = input('(1) Non-advanced solver (large errors might occur for advanced DDAEs)\n(2) Advanced solver (solution interval must not be too big)\n');
                switch choice
                    case 1
                        [t,x] = solve_lddae(E,A,B,f,tau,phi,tspan,options);
                    case 2
                        if isa(tau,'function_handle') || numel(tau)>1
                            error('The solver for advanced DDAEs only works for single constant delay.');
                        end
                        [t,x] = solve_advanced_lddae(E,A,B,f,tau,phi,tspan,options);
                end
            case 2
                F1 = data{1};
                F2 = data{2};
                tau = data{3};
                phi = data{4};
                [t,x]=solve_nddae(F1,F2,tau,phi,tspan,options);
        end
    case 2
        E = data{1};
        A = data{2};
        B = data{3};
        f = data{4};
        tau = data{5};
        phi = data{6};
        if isa(tau,'function_handle') || numel(tau)>1
            error('The solver for DDAEs with shift index bigger than zero only works for single constant delay.');
        end
        [t,x] = solve_shifted_lddae(E,A,B,f,tau,phi,tspan,options);
end