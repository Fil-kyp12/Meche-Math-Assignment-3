%using the explicit midpoint approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB,num_evals] = explicit_midpoint_step(rate_func_in,t,XA,h)

    % First evaluation: slope at the beginning of the interval
    k1 = rate_func_in(t, XA);

    % Second evaluation: slope at the midpoint
    k2 = rate_func_in(t + h/2, XA + (h/2)*k1);

    % Update X using the midpoint slope
    XB = XA + h*k2;

    % The rate function was called twice
    num_evals = 2;

end