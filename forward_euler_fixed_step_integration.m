%Runs numerical integration using forward Euler approximation
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
    forward_euler_fixed_step_integration(rate_func_in,tspan,X0,h_ref)

    t_start = tspan(1);
    t_end = tspan(2);

    N = ceil((t_end - t_start)/h_ref);

    h_avg = (t_end - t_start)/N;

    t_list = zeros(N+1,1);
    X_list = zeros(N+1,length(X0));

    t_list(1) = t_start;
    X_list(1,:) = X0(:)';

    num_evals = 0;

    for n = 1:N

        [X_next,evals] = forward_euler_step( ...
            rate_func_in,t_list(n),X_list(n,:),h_avg);

        t_list(n+1) = t_list(n) + h_avg;
        X_list(n+1,:) = X_next(:)';

        num_evals = num_evals + evals;

    end

end