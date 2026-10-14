function [XB, num_evals] = implicit_midpoint_step (rate_func_in, t, XA, h)
    G = @(X_next) XA + h*rate_func_in(t + h/2, 0.5*(XA + X_next)) - X_next;
    [XB, ~, num_evals] = multi_newton_solver(G,XA);
end