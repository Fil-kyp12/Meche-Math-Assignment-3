function [XB, num_evals] = backward_euler_step(rate_func_in, t, XA, h)
    G = @(X_next) XA + h*rate_func_in(t+h,X_next)-X_next;
    [XB, ~, num_evals] = multi_newton_solver(G,XA);
end 