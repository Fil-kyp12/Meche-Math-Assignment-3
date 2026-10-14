%% APPLIED MATHEMATICS FOR ENGINEERS
% Numerical ODE methods: all plots, error tables, and stability tests
% Requires your existing step functions and fixed_step_integration.m
% MATLAB R2016b or newer is recommended.

clearvars;
close all;
clc;
format short e;

%% 1. SETTINGS

tspan = [0, 10];
h_list = [0.3, 0.1, 0.03];

% Step sizes for error and convergence studies
h_values = logspace(log10(0.005), log10(0.5), 10);

% Methods -- uses your existing functions
methods = {@forward_euler_step, @explicit_midpoint_step, ...
           @backward_euler_step, @implicit_midpoint_step};

names = {'Forward Euler', 'Explicit Midpoint', ...
         'Backward Euler', 'Implicit Midpoint'};

% LaTeX-compatible display names
latex_names = {'Forward Euler', 'Explicit Midpoint', ...
               'Backward Euler', 'Implicit Midpoint'};

% Test problems and their exact solutions
rates = {@rate_func01, @rate_func02};
solutions = {@solution01, @solution02};
test_names = {'Test 1', 'Test 2'};

nMethods = numel(methods);
nH = numel(h_values);
nTests = numel(rates);

% Consistent plotting defaults
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultAxesLineWidth', 1.1);
set(groot, 'defaultLineLineWidth', 1.8);
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

colors = lines(nMethods);
markers = {'o', 's', '^', 'd'};

%% 2. EXAMPLE SOLUTION PLOTS WITH ERROR PANELS

tspan = [0, 10];
X0 = solution01(tspan(1));
h_list = [0.3, 0.1, 0.03];

t_true = linspace(tspan(1), tspan(2), 2000);
x_true = solution01(t_true);

% Distinct colours for the three step sizes
stepColors = [0.0000 0.4470 0.7410;
              0.8500 0.3250 0.0980;
              0.4660 0.6740 0.1880];

for m = 1:numel(methods)

    figure('Name', names{m}, 'NumberTitle', 'off', ...
           'Color', 'w');

    % Two panels: solution and absolute error
    tiledlayout(2,1, 'TileSpacing', 'compact', ...
                    'Padding', 'compact');

    %% Top panel: numerical and exact solutions
    ax1 = nexttile;
    hold(ax1, 'on');

    plot(ax1, t_true, x_true, 'k-', ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Exact: $x(t)=\cos(t)$');

    for k = 1:numel(h_list)

        [t_list, X_list, h_avg] = fixed_step_integration( ...
            @rate_func01, methods{m}, tspan, X0, h_list(k));

        % Draw the complete numerical curve
        plot(ax1, t_list, X_list, '-', ...
            'Color', stepColors(k,:), ...
            'LineWidth', 1.5, ...
            'HandleVisibility', 'off');

        % Sparse markers so small h does not clutter the graph
        markerSpacing = max(1, round(0.3/h_avg));
        idx = 1:markerSpacing:numel(t_list);

        plot(ax1, t_list(idx), X_list(idx), 'o', ...
            'Color', stepColors(k,:), ...
            'MarkerFaceColor', stepColors(k,:), ...
            'MarkerSize', 4, ...
            'LineStyle', 'none', ...
            'DisplayName', sprintf('$h=%.4f$', h_avg));
    end

    ylabel('$x(t)$', 'Interpreter', 'latex');
    xlabel('$t$', 'Interpreter', 'latex');
    title([names{m} ': numerical versus exact solution'], ...
          'Interpreter', 'latex');
    legend(ax1, 'Location', 'best', 'Box', 'off');
    grid(ax1, 'on');
    box(ax1, 'on');
    xlim(ax1, tspan);

    %% Bottom panel: absolute error
    ax2 = nexttile;
    hold(ax2, 'on');

    for k = 1:numel(h_list)

        [t_list, X_list, h_avg] = fixed_step_integration( ...
            @rate_func01, methods{m}, tspan, X0, h_list(k));

        x_exact_at_nodes = solution01(t_list);
        abs_error = abs(X_list(:) - x_exact_at_nodes(:));

        plot(ax2, t_list, abs_error, '-', ...
            'Color', stepColors(k,:), ...
            'LineWidth', 1.7, ...
            'DisplayName', sprintf('$h=%.4f$', h_avg));
    end

    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$|x_{\mathrm{num}}-x_{\mathrm{exact}}|$', ...
           'Interpreter', 'latex');
    title('Absolute error', 'Interpreter', 'latex');
    legend(ax2, 'Location', 'best', 'Box', 'off');
    grid(ax2, 'on');
    box(ax2, 'on');
    xlim(ax2, tspan);

    % Shared overall title
    sgtitle([names{m} ' ($x''=-5x+5\cos(t)-\sin(t)$)'], ...
            'Interpreter', 'latex');
end

%% 3. LOCAL TRUNCATION ERRORS
% Start each one-step calculation from the exact solution at t_ref.
% Local error = || one numerical step - exact value at t_ref+h ||_2.

t_ref = 0.5;
localErr = zeros(nTests, nMethods, nH);

for q = 1:nTests
    for m = 1:nMethods
        for k = 1:nH
            h = h_values(k);
            Xref = solutions{q}(t_ref);

            [Xnum, ~] = methods{m}(rates{q}, t_ref, Xref, h);
            Xexact = solutions{q}(t_ref + h);

            localErr(q,m,k) = norm(Xnum(:) - Xexact(:), 2);
        end
    end
end

% Estimate convergence exponents from log(error) = p*log(h) + c.
pLocal = nan(nTests, nMethods);

for q = 1:nTests
    for m = 1:nMethods
        pLocal(q,m) = estimate_slope( ...
            h_values, squeeze(localErr(q,m,:))');
    end
end


LocalErrorTable = table( ...
    names', pLocal(1,:)', pLocal(2,:)', ...
    'VariableNames', ...
    {'Method', 'Test1_p_minus', 'Test2_p_minus'});

disp(' ');
disp('Table 1. Local Truncation Error Scaling vs. Step Size');
disp('Convergence order p (-), dimensionless');
disp(LocalErrorTable);


%% 4. LOCAL ERROR PLOTS: ALL METHODS, TEST 1

figure('Name', 'Local error - all methods', 'NumberTitle', 'off');
hold on;

for m = 1:nMethods
    loglog(h_values, squeeze(localErr(1,m,:)), ...
        [markers{m} '-'], ...
        'Color', colors(m,:), ...
        'MarkerSize', 6, ...
        'DisplayName', latex_names{m});
end

xlabel('Step size, $h$');
ylabel('Local error, $\|X_{n+1}^{\mathrm{num}}-X(t+h)\|_2$');
title('Local truncation error: Test 1');
legend('Location', 'best', 'Box', 'off');
grid on; box on;
hold off;


%% 5. LOCAL ERROR: EXPLICIT METHODS AND REFERENCE INCREMENT

figure('Name', 'Local error - explicit methods', ...
       'NumberTitle', 'off', ...
       'Color', 'w');

ax = axes;
hold(ax, 'on');

explicit_idx = [1, 2];

for j = 1:numel(explicit_idx)
    m = explicit_idx(j);
    err = squeeze(localErr(1,m,:))';
    
    % Measured local errors
    loglog(ax, h_values, max(err, realmin), ...
        [markers{m} '-'], ...
        'Color', colors(m,:), ...
        'MarkerSize', 7, ...
        'LineWidth', 1.8, ...
        'DisplayName', latex_names{m});

    % Fitted convergence line using positive, finite errors
    valid = isfinite(err) & err > 0 & ...
            isfinite(h_values) & h_values > 0;

    if nnz(valid) >= 2
        p = estimate_slope(h_values(valid), err(valid));

        % Fit only over the valid data range
        h_fit = h_values(valid);
        err_fit = err(valid);

        fit_y = fitted_values(h_fit, err_fit, p);

        loglog(ax, h_fit, fit_y, '--', ...
            'Color', colors(m,:), ...
            'LineWidth', 2, ...
            'HandleVisibility', 'off');
    end
end

% Exact-solution increment for comparison
exact_increment = abs(solution01(t_ref + h_values) ...
                     - solution01(t_ref));

loglog(ax, h_values, max(exact_increment, realmin), 'k:', ...
    'LineWidth', 2.2, ...
    'DisplayName', '$|x(t+h)-x(t)|$');

% Labels and title
xlabel(ax, 'Step size, $h$', ...
    'Interpreter', 'latex', 'FontSize', 14);

ylabel(ax, 'Local truncation error', ...
    'Interpreter', 'latex', 'FontSize', 14);

title(ax, 'Local Truncation Error Scaling: Explicit Methods', ...
    'Interpreter', 'latex', 'FontSize', 15);

% LaTeX tick labels and readable axes
set(ax, 'TickLabelInterpreter', 'latex', ...
    'FontSize', 12, ...
    'LineWidth', 1.1);

legend(ax, 'Location', 'best', ...
    'Interpreter', 'latex', ...
    'FontSize', 11, ...
    'Box', 'off');

grid(ax, 'on');
box(ax, 'on');
hold(ax, 'off');


%% 6. GLOBAL ENDPOINT ERRORS AND FUNCTION EVALUATIONS
% Global error is measured at t = tspan(2).
% Calls are the rate-function evaluations reported by your integrator.

globalErr = zeros(nTests, nMethods, nH);
numCalls = zeros(nTests, nMethods, nH);
actualH = zeros(nTests, nMethods, nH);

for q = 1:nTests
    X0q = solutions{q}(tspan(1));
    Xend_exact = solutions{q}(tspan(2));

    for m = 1:nMethods
        for k = 1:nH
            [t_list, X_list, h_avg, evals] = ...
                fixed_step_integration( ...
                    rates{q}, methods{m}, tspan, X0q, h_values(k));

            Xend_num = X_list(end,:)';

            globalErr(q,m,k) = norm( ...
                Xend_num(:) - Xend_exact(:), 2);

            numCalls(q,m,k) = evals;
            actualH(q,m,k) = h_avg;
        end
    end
end

%% 7. GLOBAL ERROR CONVERGENCE ORDERS

pGlobal = nan(nTests, nMethods);
pCalls = nan(nTests, nMethods);

for q = 1:nTests
    for m = 1:nMethods
        err = squeeze(globalErr(q,m,:))';
        hs = squeeze(actualH(q,m,:))';
        calls = squeeze(numCalls(q,m,:))';

        pGlobal(q,m) = estimate_slope(hs, err);

        % Error versus calls generally scales like calls^(-p).
        pCalls(q,m) = -estimate_slope(calls, err);
    end
end

GlobalStepTable = table( ...
    names', pGlobal(1,:)', pGlobal(2,:)', ...
    'VariableNames', {'Method', 'Test1_p_minus', 'Test2_p_minus'});

GlobalCallsTable = table( ...
    names', pCalls(1,:)', pCalls(2,:)', ...
    'VariableNames', {'Method', 'Test1_p_minus', 'Test2_p_minus'});

disp(' ');
disp('===== GLOBAL ERROR ORDERS: ERROR VS STEP SIZE =====');
disp(GlobalStepTable);

disp(' ');
disp('===== GLOBAL ERROR ORDERS: ERROR VS FUNCTION CALLS =====');
disp(GlobalCallsTable);

%% 8. GLOBAL ERROR VS STEP SIZE: ALL METHODS

figure('Name', 'Global error versus step size', ...
       'NumberTitle', 'off');
hold on;

for m = 1:nMethods
    err = squeeze(globalErr(1,m,:))';
    hs = squeeze(actualH(1,m,:))';

    loglog(hs, max(err, realmin), [markers{m} '-'], ...
        'Color', colors(m,:), ...
        'MarkerSize', 6, ...
        'DisplayName', latex_names{m});
end

xlabel('Step size, $h$');
ylabel('Global endpoint error, $\|X_N-X(T)\|_2$');
title('Global endpoint error: Test 1');
legend('Location', 'best', 'Box', 'off');
grid on; box on;
hold off;

%% 9. GLOBAL ERROR VS STEP SIZE: EXPLICIT METHODS WITH FITS

figure('Name', 'Global error - explicit methods', ...
       'NumberTitle', 'off');
hold on;

for m = explicit_idx
    err = squeeze(globalErr(1,m,:))';
    hs = squeeze(actualH(1,m,:))';

    loglog(hs, max(err, realmin), [markers{m} '-'], ...
        'Color', colors(m,:), ...
        'MarkerSize', 6, ...
        'DisplayName', latex_names{m});

    valid = err > 0;
    if nnz(valid) >= 2
        p = estimate_slope(hs(valid), err(valid));
        loglog(hs, fitted_values(hs, err, p), '--', ...
            'Color', colors(m,:), ...
            'HandleVisibility', 'off');
    end
end

xlabel('Step size, $h$');
ylabel('Global endpoint error, $\|X_N-X(T)\|_2$');
title('Global error scaling: explicit methods');
legend('Location', 'best', 'Box', 'off');
grid on; box on;
hold off;

%% 10. GLOBAL ERROR VS RATE-FUNCTION CALLS: ALL METHODS

figure('Name', 'Global error versus function calls', ...
       'NumberTitle', 'off');
hold on;

for m = 1:nMethods
    err = squeeze(globalErr(1,m,:))';
    calls = squeeze(numCalls(1,m,:))';

    loglog(calls, max(err, realmin), [markers{m} '-'], ...
        'Color', colors(m,:), ...
        'MarkerSize', 6, ...
        'DisplayName', latex_names{m});
end

xlabel('Rate-function evaluations, $N_f$');
ylabel('Global endpoint error, $\|X_N-X(T)\|_2$');
title('Global error versus computational work: Test 1');
legend('Location', 'best', 'Box', 'off');
grid on; box on;
hold off;

%% 11. GLOBAL ERROR VS FUNCTION CALLS: EXPLICIT METHODS WITH FITS

figure('Name', 'Global error versus calls - explicit methods', ...
       'NumberTitle', 'off');
hold on;

for m = explicit_idx
    err = squeeze(globalErr(1,m,:))';
    calls = squeeze(numCalls(1,m,:))';

    loglog(calls, max(err, realmin), [markers{m} '-'], ...
        'Color', colors(m,:), ...
        'MarkerSize', 6, ...
        'DisplayName', latex_names{m});

    valid = err > 0;
    if nnz(valid) >= 2
        p = estimate_slope(calls(valid), err(valid));
        fit_y = fitted_values(calls, err, p);
        loglog(calls, fit_y, '--', ...
            'Color', colors(m,:), ...
            'HandleVisibility', 'off');
    end
end

xlabel('Rate-function evaluations, $N_f$');
ylabel('Global endpoint error, $\|X_N-X(T)\|_2$');
title('Error versus computational work: explicit methods');
legend('Location', 'best', 'Box', 'off');
grid on; box on;
hold off;

%% 12. INTEGRATION STABILITY TEST: h_ref = 0.38

stability_h = 0.38;
stability_tspan = [0, 20];
X0_stability = solution01(stability_tspan(1));

t_stability = linspace(stability_tspan(1), ...
                       stability_tspan(2), 2000);
x_stability = solution01(t_stability);

% One figure containing four vertically stacked subplots
figure('Name', 'Integration Stability, h = 0.38', ...
       'NumberTitle', 'off', ...
       'Color', 'w');

tiledlayout(4,1, 'TileSpacing', 'compact', ...
                 'Padding', 'compact');

for m = 1:nMethods

    ax = nexttile;
    hold(ax, 'on');

    % Analytical solution
    plot(ax, t_stability, x_stability, 'k-', ...
        'LineWidth', 2.2, ...
        'DisplayName', '$x(t)=\cos(t)$');

    % Numerical solution using h = 0.38
    [t_num, X_num, h_avg] = fixed_step_integration( ...
        @rate_func01, methods{m}, stability_tspan, ...
        X0_stability, stability_h);

    plot(ax, t_num, X_num, 'o-', ...
        'Color', colors(m,:), ...
        'MarkerSize', 3, ...
        'LineWidth', 1.3, ...
        'DisplayName', sprintf('Numerical, $h=%.3f$', h_avg));

    xlabel(ax, '$t$', 'Interpreter', 'latex');
    ylabel(ax, '$x(t)$', 'Interpreter', 'latex');

    title(ax, sprintf('%s: $h_{\\mathrm{ref}}=%.2f$', ...
          names{m}, stability_h), ...
          'Interpreter', 'latex');

    legend(ax, 'Location', 'best', 'Box', 'off');

    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, stability_tspan);

end

sgtitle('Integration Stability Comparison', ...
        'Interpreter', 'latex');
%% 13. INTEGRATION STABILITY TEST: h_ref = 0.45

stability_h = 0.45;
stability_tspan = [0, 20];
X0_stability = solution01(stability_tspan(1));

t_stability = linspace(stability_tspan(1), ...
                       stability_tspan(2), 2000);
x_stability = solution01(t_stability);

% One figure containing four vertically stacked subplots
figure('Name', 'Integration Stability, h = 0.45', ...
       'NumberTitle', 'off', ...
       'Color', 'w');

tiledlayout(4,1, 'TileSpacing', 'compact', ...
                 'Padding', 'compact');

for m = 1:nMethods

    ax = nexttile;
    hold(ax, 'on');

    % Analytical solution
    plot(ax, t_stability, x_stability, 'k-', ...
        'LineWidth', 2.2, ...
        'DisplayName', '$x(t)=\cos(t)$');

    % Numerical solution using h = 0.45
    [t_num, X_num, h_avg] = fixed_step_integration( ...
        @rate_func01, methods{m}, stability_tspan, ...
        X0_stability, stability_h);

    plot(ax, t_num, X_num, 'o-', ...
        'Color', colors(m,:), ...
        'MarkerSize', 3, ...
        'LineWidth', 1.3, ...
        'DisplayName', sprintf('Numerical, $h=%.3f$', h_avg));

    xlabel(ax, '$t$', 'Interpreter', 'latex');
    ylabel(ax, '$x(t)$', 'Interpreter', 'latex');

    title(ax, sprintf('%s: $h_{\\mathrm{ref}}=%.2f$', ...
          names{m}, stability_h), ...
          'Interpreter', 'latex');

    legend(ax, 'Location', 'best', 'Box', 'off');

    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, stability_tspan);

end

sgtitle('Integration Stability Comparison, $h_{\\mathrm{ref}}=0.45$', ...
        'Interpreter', 'latex');

%% LOCAL HELPER FUNCTIONS
% These helpers do not replace or modify any of your ODE methods.

function p = estimate_slope(x, y)
    x = x(:);
    y = y(:);

    valid = isfinite(x) & isfinite(y) & x > 0 & y > 0;
    x = x(valid);
    y = y(valid);

    if numel(x) < 2
        p = NaN;
        return;
    end

    coeff = polyfit(log(x), log(y), 1);
    p = coeff(1);
end

function yfit = fitted_values(x, y, p)
    x = x(:)';
    y = y(:)';

    valid = isfinite(x) & isfinite(y) & x > 0 & y > 0;

    if nnz(valid) < 2 || ~isfinite(p)
        yfit = max(y, realmin);
        return;
    end

    % Fit a line with the estimated slope p in log-log coordinates.
    intercept = mean(log(y(valid)) - p .* log(x(valid)));
    yfit = exp(intercept) .* x.^p;
    yfit = max(yfit, realmin);
end