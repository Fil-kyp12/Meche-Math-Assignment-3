methods = {@forward_euler_step, @explicit_midpoint_step, @backward_euler_step, @implicit_midpoint_step};
for m = 1:4
    [t, X] = fixed_step_integration(@rate_func01, methods{m}, [0 10], 1, 0.1);
    figure; plot(t, X, t, cos(t), '--'); title(func2str(methods{m}))
end