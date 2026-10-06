%Runs numerical integration using explicit midpoint approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg,num_evals] = ...
    explicit_midpoint_fixed_step_integration(rate_func_in,tspan,X0,h_ref)

    % Get integration endpoints
    t_start = tspan(1);
    t_end = tspan(2);

    % Find smallest N such that h < h_ref
    N = ceil((t_end - t_start) / h_ref);

    % Actual step size
    h = (t_end - t_start) / N;

    % Average step size
    h_avg = h;

    % Initialize time and solution arrays
    t_list = zeros(N+1,1);
    X_list = zeros(N+1,length(X0));

    % Initial condition
    t_list(1) = t_start;
    X_list(1,:) = X0(:)';

    % Initialize number of function evaluations
    num_evals = 0;

    % Perform the integration
    for n = 1:N

        % Current time and X
        t = t_list(n);
        XA = X_list(n,:)';

        % Take one explicit midpoint step
        [XB, evals] = explicit_midpoint_step(rate_func_in,t,XA,h);

        % Store result
        t_list(n+1) = t + h;
        X_list(n+1,:) = XB(:)';

        % Update evaluation count
        num_evals = num_evals + evals;

    end

end