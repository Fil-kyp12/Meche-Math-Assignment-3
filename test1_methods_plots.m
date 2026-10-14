tspan  = [0, 10];
X0     = solution01(tspan(1));
h_list = [0.3, 0.1, 0.03];  
methods = {@forward_euler_step, @explicit_midpoint_step, ...
    @backward_euler_step, @implicit_midpoint_step};
names   = {'Forward Euler', 'Explicit Midpoint', ...
    'Backward Euler', 'Implicit Midpoint'};
% Analytical solution on a fine grid
t_true = linspace(tspan(1), tspan(2), 1000);
x_true = solution01(t_true);
for m = 1:length(methods)
figure; hold on
plot(t_true, x_true, 'k', 'LineWidth', 2, 'DisplayName', 'Analytical: cos(t)')

for k = 1:length(h_list)
    [t_list, X_list, h_avg] = fixed_step_integration( ...
        @rate_func01, methods{m}, tspan, X0, h_list(k));
    plot(t_list, X_list, 'o-', 'MarkerSize', 3, ...
        'DisplayName', sprintf('h = %.3f', h_avg))
end

xlabel('t'); ylabel('x(t)')
title([names{m} ': Numerical vs. Analytical Solution (Test 1)'])
legend('Location', 'best'); grid on
hold off
end