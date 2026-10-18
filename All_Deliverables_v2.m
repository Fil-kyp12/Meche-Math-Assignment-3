%% APPLIED MATHEMATICS FOR ENGINEERS
% Numerical ODE methods: ALL 15 deliverables (plots, error tables, stability tests)
% Same layout as All_Deliverables.m, corrected to match the handout + Canvas rules.
% Requires your existing step functions, fixed_step_integration.m,
% multi_newton_solver.m, approximate_jacobian.m, rate_func01/02.m, solution01/02.m
% (all in the same folder; none of them are modified).
% MATLAB R2020a or newer (exportgraphics). Full run: roughly 5-15 minutes.
% Output goes to the folder  deliverables_out  (PNG plots at 600 dpi + table files).
%
% Deliverable numbers (Canvas):
%  #1-4  example plots     #5-6  stability     #7-8  local error plots   #9  local table
%  #10-11 global vs h      #12   table         #13-14 global vs calls    #15 table

clearvars;
close all;
clc;
format short e;

FAST = false;     % true = quick test run (smaller step-size range, ~1 minute)

%% 1. SETTINGS

tspan = [0, 10];
h_list = [0.35, 0.2, 0.1];           % example plots: all well below the 0.38 stability limit

% Step sizes for the error studies (handout: local ~100 points, 1e-5 .. 1e1)
h_local = logspace(-5, 1, 100);
if FAST
    h_values = logspace(-3, 0, 13);
else
    h_values = logspace(-4, 0, 25);  % global study, four decades
end

% Fit ONLY the straight part of each curve (above round-off, below instability)
fit_lo = 1e-4;  fit_hi = 1e-2;
if FAST, fit_lo = 1e-3; end

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
nHL = numel(h_local);
nTests = numel(rates);

explicit_idx = [1, 2];               % Forward Euler, Explicit Midpoint

% Consistent plotting defaults: LaTeX for text, title, legend AND tick labels
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

colors = [0.84 0.15 0.16;    % Forward Euler
          0.12 0.47 0.71;    % Explicit Midpoint
          0.58 0.40 0.74;    % Backward Euler
          0.17 0.63 0.17];   % Implicit Midpoint
markers = {'o', 's', 'd', '^'};
% FE/BE and EM/IM have almost identical error curves. Explicit methods (1,2) get big markers and are
% drawn first; implicit methods (3,4) get small markers and are drawn on top, so both stay visible.
msize_nofit = [8 8 4.5 4.5];

outdir = 'deliverables_out';
if ~exist(outdir, 'dir'), mkdir(outdir); end

fit_title = 'Fit window: $10^{-4}\le h\le10^{-2}$';

%% 2. EXAMPLE SOLUTION PLOTS (#1-#4): one plot per method

X0 = solution01(tspan(1));

t_true = linspace(tspan(1), tspan(2), 2000);
x_true = solution01(t_true);

% Distinct colours / styles for the three step sizes
stepColors = [0.85 0.37 0.01;
              0.11 0.62 0.47;
              0.46 0.44 0.70];
stepMarkers = {'o', 's', 'none'};
stepStyles  = {'-', '-', '--'};

example_names = {'01_forward_euler_example', '02_explicit_midpoint_example', ...
                 '03_backward_euler_example', '04_implicit_midpoint_example'};

for m = 1:numel(methods)

    f = figure('Name', names{m}, 'NumberTitle', 'off', ...
               'Color', 'w', 'Position', [100 100 900 650]);
    ax1 = axes;
    hold(ax1, 'on');

    plot(ax1, t_true, x_true, 'k-', 'LineWidth', 3, ...
        'DisplayName', 'Analytical: $x(t)=\cos(t)$');

    for k = 1:numel(h_list)

        [t_list, X_list, h_avg] = fixed_step_integration( ...
            @rate_func01, methods{m}, tspan, X0, h_list(k));

        plot(ax1, t_list, X_list(:,1), ...
            'LineStyle', stepStyles{k}, ...
            'Color', stepColors(k,:), ...
            'Marker', stepMarkers{k}, ...
            'MarkerSize', 7, ...
            'MarkerFaceColor', stepColors(k,:), ...
            'MarkerEdgeColor', stepColors(k,:), ...
            'LineWidth', 2, ...
            'DisplayName', sprintf('Numerical, $h=%.3f$', h_avg));
    end

    xlabel(ax1, 'Time, $t$ (-)', 'FontSize', 18);
    ylabel(ax1, 'Solution, $x(t)$ (-)', 'FontSize', 18);
    title(ax1, [names{m} ': numerical vs.\ analytical solution'], 'FontSize', 18);
    legend(ax1, 'Location', 'northeast', 'FontSize', 15);
    grid(ax1, 'on');
    box(ax1, 'on');
    xlim(ax1, tspan);
    ylim(ax1, [-1.15 2.6]);              % head-room so the legend covers no data

    save_figure(f, fullfile(outdir, example_names{m}));
end

%% 3. LOCAL TRUNCATION ERRORS
% Start each one-step calculation from the exact solution at t_ref.
% Local error = || one numerical step - exact value at t_ref+h ||_2.

t_ref = 0.492;                           % nonzero, x'(t_ref) is nonzero
localErr = zeros(nTests, nMethods, nHL);
diffX = zeros(nTests, nHL);              % |X(t+h) - X(t)| reference

for q = 1:nTests
    Xref = solutions{q}(t_ref);
    for k = 1:nHL
        diffX(q,k) = norm(solutions{q}(t_ref + h_local(k)) - Xref, 2);
    end
    for m = 1:nMethods
        for k = 1:nHL
            h = h_local(k);

            [Xnum, ~] = methods{m}(rates{q}, t_ref, Xref, h);
            Xexact = solutions{q}(t_ref + h);

            localErr(q,m,k) = norm(Xnum(:) - Xexact(:), 2);
        end
    end
end

% Estimate convergence exponents from log(error) = p*log(h) + c  (fit window only)
win_local = h_local >= fit_lo*(1-1e-9) & h_local <= fit_hi*(1+1e-9);
pLocal = nan(nTests, nMethods);
pDiff = nan(1, nTests);

for q = 1:nTests
    pDiff(q) = estimate_slope(h_local(win_local), diffX(q,win_local));
    for m = 1:nMethods
        pLocal(q,m) = estimate_slope( ...
            h_local(win_local), squeeze(localErr(q,m,win_local))');
    end
end

LocalErrorTable = table( ...
    names', pLocal(1,:)', pLocal(2,:)', ...
    'VariableNames', ...
    {'Method', 'Test1_p_minus', 'Test2_p_minus'});

disp(' ');
disp('Table 1 (#9). Local Truncation Error Scaling vs. Step Size');
disp('Convergence order p (-), dimensionless');
disp(LocalErrorTable);
fprintf('|X(t+h)-X(t)| row: Test1 p = %.3f, Test2 p = %.3f\n', pDiff(1), pDiff(2));

%% 4. LOCAL ERROR PLOT: ALL METHODS, NO FITS (#8)

f = figure('Name', 'Local error - all methods', 'NumberTitle', 'off', ...
           'Color', 'w', 'Position', [100 100 900 680]);
hold on;

for m = 1:nMethods
    loglog(h_local, squeeze(localErr(1,m,:)), ...
        ['-' markers{m}], ...
        'Color', colors(m,:), ...
        'MarkerSize', msize_nofit(m), ...
        'MarkerFaceColor', colors(m,:), ...
        'MarkerEdgeColor', colors(m,:), ...
        'LineWidth', 1.4, ...
        'DisplayName', latex_names{m});
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Step size, $h$ (-)', 'FontSize', 18);
ylabel('Local truncation error, $\epsilon_{\mathrm{loc}}$ (-)', 'FontSize', 18);
title('Local truncation error (all methods, Test 1)', 'FontSize', 18);
legend('Location', 'best', 'FontSize', 15);
grid on; box on;
hold off;

save_figure(f, fullfile(outdir, '08_local_error_loglog_all_methods'));

%% 5. LOCAL ERROR: EXPLICIT METHODS, FIT LINES AND REFERENCE INCREMENT (#7)

f = figure('Name', 'Local error - explicit methods', 'NumberTitle', 'off', ...
           'Color', 'w', 'Position', [100 100 900 680]);

ax = axes;
hold(ax, 'on');

for j = 1:numel(explicit_idx)
    m = explicit_idx(j);
    err = squeeze(localErr(1,m,:))';

    % Measured local errors (data points only)
    loglog(ax, h_local, err, markers{m}, ...
        'Color', colors(m,:), ...
        'MarkerSize', 4, ...
        'MarkerFaceColor', colors(m,:), ...
        'MarkerEdgeColor', colors(m,:), ...
        'LineStyle', 'none', ...
        'DisplayName', latex_names{m});

    % Fit line over the straight part only
    [h_fit, fit_y, p] = fit_line(h_local, err, win_local);
    loglog(ax, h_fit, fit_y, '-', 'Color', 'w', 'LineWidth', 6, 'HandleVisibility', 'off');  % halo
    loglog(ax, h_fit, fit_y, '--', ...
        'Color', colors(m,:), ...
        'LineWidth', 2.5, ...
        'DisplayName', sprintf('%s fit, $p=%.2f$', latex_names{m}, p));
end

% Exact-solution increment |X(t+h)-X(t)| with its own fit line
loglog(ax, h_local, diffX(1,:), 'x', ...
    'Color', 'k', 'MarkerSize', 5, 'LineStyle', 'none', ...
    'DisplayName', '$|X(t+h)-X(t)|$');
[h_fit, fit_y, p] = fit_line(h_local, diffX(1,:), win_local);
loglog(ax, h_fit, fit_y, '-', 'Color', 'w', 'LineWidth', 6, 'HandleVisibility', 'off');  % halo
loglog(ax, h_fit, fit_y, 'k--', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('$|X(t+h)-X(t)|$ fit, $p=%.2f$', p));

set(ax, 'XScale', 'log', 'YScale', 'log');
xlabel(ax, 'Step size, $h$ (-)', 'FontSize', 18);
ylabel(ax, 'Local truncation error, $\epsilon_{\mathrm{loc}}$ (-)', 'FontSize', 18);
title(ax, 'Local truncation error (explicit methods, Test 1)', 'FontSize', 18);
lg = legend(ax, 'Location', 'best', 'FontSize', 13);
title(lg, fit_title);
grid(ax, 'on');
box(ax, 'on');
hold(ax, 'off');

save_figure(f, fullfile(outdir, '07_local_error_loglog_fit_explicit'));

%% 6. GLOBAL ENDPOINT ERRORS AND FUNCTION EVALUATIONS
% Global error is measured at t = tspan(2).
% Calls are the rate-function evaluations reported by your integrator.

globalErr = zeros(nTests, nMethods, nH);
numCalls = zeros(nTests, nMethods, nH);
actualH = zeros(nTests, nMethods, nH);

for q = 1:nTests
    X0q = solutions{q}(tspan(1));

    for m = 1:nMethods
        fprintf('Global sweep: %s, %s ...\n', test_names{q}, names{m});
        for k = 1:nH
            [t_list, X_list, h_avg, evals] = ...
                fixed_step_integration( ...
                    rates{q}, methods{m}, tspan, X0q, h_values(k));

            Xend_exact = solutions{q}(t_list(end));
            Xend_num = X_list(end,:)';

            globalErr(q,m,k) = norm( ...
                Xend_num(:) - Xend_exact(:), 2);

            numCalls(q,m,k) = evals;
            actualH(q,m,k) = h_avg;
        end
    end
end

%% 7. GLOBAL ERROR CONVERGENCE ORDERS (fit window only)

win_glob = h_values >= fit_lo*(1-1e-9) & h_values <= fit_hi*(1+1e-9);
pGlobal = nan(nTests, nMethods);
pCalls = nan(nTests, nMethods);

for q = 1:nTests
    for m = 1:nMethods
        err = squeeze(globalErr(q,m,:))';
        hs = squeeze(actualH(q,m,:))';
        calls = squeeze(numCalls(q,m,:))';

        pGlobal(q,m) = estimate_slope(hs(win_glob), err(win_glob));

        % Error vs calls scales like calls^p (p is negative: -1 for Euler, -2 for midpoint)
        pCalls(q,m) = estimate_slope(calls(win_glob), err(win_glob));
    end
end

GlobalStepTable = table( ...
    names', pGlobal(1,:)', pGlobal(2,:)', ...
    'VariableNames', {'Method', 'Test1_p_minus', 'Test2_p_minus'});

GlobalCallsTable = table( ...
    names', pCalls(1,:)', pCalls(2,:)', ...
    'VariableNames', {'Method', 'Test1_p_minus', 'Test2_p_minus'});

disp(' ');
disp('===== (#12) GLOBAL ERROR ORDERS: ERROR VS STEP SIZE =====');
disp(GlobalStepTable);

disp(' ');
disp('===== (#15) GLOBAL ERROR ORDERS: ERROR VS FUNCTION CALLS =====');
disp(GlobalCallsTable);

%% 8. GLOBAL ERROR VS STEP SIZE: ALL METHODS, NO FITS (#11)

f = figure('Name', 'Global error versus step size', 'NumberTitle', 'off', ...
           'Color', 'w', 'Position', [100 100 900 680]);
hold on;

for m = 1:nMethods
    err = squeeze(globalErr(1,m,:))';
    hs = squeeze(actualH(1,m,:))';

    loglog(hs, err, ['-' markers{m}], ...
        'Color', colors(m,:), ...
        'MarkerSize', msize_nofit(m), ...
        'MarkerFaceColor', colors(m,:), ...
        'MarkerEdgeColor', colors(m,:), ...
        'LineWidth', 1.4, ...
        'DisplayName', latex_names{m});
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Step size, $h$ (-)', 'FontSize', 18);
ylabel('Global truncation error, $\epsilon_{\mathrm{glob}}$ (-)', 'FontSize', 18);
title('Global truncation error vs.\ step size (all methods, Test 1)', 'FontSize', 18);
legend('Location', 'best', 'FontSize', 15);
grid on; box on;
hold off;

save_figure(f, fullfile(outdir, '11_global_error_vs_h_all_methods'));

%% 9. GLOBAL ERROR VS STEP SIZE: EXPLICIT METHODS WITH FITS (#10)

f = figure('Name', 'Global error - explicit methods', 'NumberTitle', 'off', ...
           'Color', 'w', 'Position', [100 100 900 680]);
hold on;

for m = explicit_idx
    err = squeeze(globalErr(1,m,:))';
    hs = squeeze(actualH(1,m,:))';

    loglog(hs, err, markers{m}, ...
        'Color', colors(m,:), ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', colors(m,:), ...
        'MarkerEdgeColor', colors(m,:), ...
        'LineStyle', 'none', ...
        'DisplayName', latex_names{m});

    [h_fit, fit_y, p] = fit_line(hs, err, win_glob);
    loglog(h_fit, fit_y, '--', ...
        'Color', colors(m,:), ...
        'LineWidth', 3, ...
        'DisplayName', sprintf('%s fit, $p=%.2f$', latex_names{m}, p));
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Step size, $h$ (-)', 'FontSize', 18);
ylabel('Global truncation error, $\epsilon_{\mathrm{glob}}$ (-)', 'FontSize', 18);
title('Global truncation error vs.\ step size (explicit methods, Test 1)', 'FontSize', 18);
lg = legend('Location', 'best', 'FontSize', 13);
title(lg, fit_title);
grid on; box on;
hold off;

save_figure(f, fullfile(outdir, '10_global_error_vs_h_fit_explicit'));

%% 10. GLOBAL ERROR VS RATE-FUNCTION CALLS: ALL METHODS, NO FITS (#14)

f = figure('Name', 'Global error versus function calls', 'NumberTitle', 'off', ...
           'Color', 'w', 'Position', [100 100 900 680]);
hold on;

for m = 1:nMethods
    err = squeeze(globalErr(1,m,:))';
    calls = squeeze(numCalls(1,m,:))';

    loglog(calls, err, ['-' markers{m}], ...
        'Color', colors(m,:), ...
        'MarkerSize', msize_nofit(m), ...
        'MarkerFaceColor', colors(m,:), ...
        'MarkerEdgeColor', colors(m,:), ...
        'LineWidth', 1.4, ...
        'DisplayName', latex_names{m});
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Rate function calls, $N_f$ (-)', 'FontSize', 18);
ylabel('Global truncation error, $\epsilon_{\mathrm{glob}}$ (-)', 'FontSize', 18);
title('Global error vs.\ rate function calls (all methods, Test 1)', 'FontSize', 18);
legend('Location', 'best', 'FontSize', 15);
grid on; box on;
hold off;

save_figure(f, fullfile(outdir, '14_global_error_vs_calls_all_methods'));

%% 11. GLOBAL ERROR VS FUNCTION CALLS: EXPLICIT METHODS WITH FITS (#13)

f = figure('Name', 'Global error versus calls - explicit methods', ...
           'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 900 680]);
hold on;

for m = explicit_idx
    err = squeeze(globalErr(1,m,:))';
    calls = squeeze(numCalls(1,m,:))';

    loglog(calls, err, markers{m}, ...
        'Color', colors(m,:), ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', colors(m,:), ...
        'MarkerEdgeColor', colors(m,:), ...
        'LineStyle', 'none', ...
        'DisplayName', latex_names{m});

    [c_fit, fit_y, p] = fit_line(calls, err, win_glob);
    loglog(c_fit, fit_y, '--', ...
        'Color', colors(m,:), ...
        'LineWidth', 3, ...
        'DisplayName', sprintf('%s fit, $p=%.2f$', latex_names{m}, p));
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Rate function calls, $N_f$ (-)', 'FontSize', 18);
ylabel('Global truncation error, $\epsilon_{\mathrm{glob}}$ (-)', 'FontSize', 18);
title('Global error vs.\ rate function calls (explicit methods, Test 1)', 'FontSize', 18);
lg = legend('Location', 'best', 'FontSize', 13);
title(lg, fit_title);
grid on; box on;
hold off;

save_figure(f, fullfile(outdir, '13_global_error_vs_calls_fit_explicit'));

%% 12. INTEGRATION STABILITY TEST: h_ref = 0.38 (#5) and h_ref = 0.45 (#6)

stability_tspan = [0, 20];
X0_stability = solution01(stability_tspan(1));

t_stability = linspace(stability_tspan(1), ...
                       stability_tspan(2), 4000);
x_stability = solution01(t_stability);

stability_h_list = [0.38, 0.45];
stability_names = {'05_stability_href_0.38', '06_stability_href_0.45'};
stability_order = [1, 3, 2, 4];     % Forward Euler, Backward Euler, Explicit Mid., Implicit Mid.

for s = 1:numel(stability_h_list)

    stability_h = stability_h_list(s);

    % One figure containing four vertically stacked subplots
    f = figure('Name', sprintf('Integration Stability, h = %.2f', stability_h), ...
               'NumberTitle', 'off', 'Color', 'w', ...
               'Position', [100 20 850 1350]);

    tiledlayout(4,1, 'TileSpacing', 'compact', ...
                     'Padding', 'compact');

    for j = 1:numel(stability_order)

        m = stability_order(j);

        ax = nexttile;
        hold(ax, 'on');

        % Analytical solution
        plot(ax, t_stability, x_stability, 'k-', ...
            'LineWidth', 3, ...
            'DisplayName', 'Analytical: $\cos(t)$');

        % Numerical solution
        [t_num, X_num, h_avg] = fixed_step_integration( ...
            @rate_func01, methods{m}, stability_tspan, ...
            X0_stability, stability_h);

        label = sprintf('Numerical, $h=%.3f$', h_avg);
        mx = max(abs(X_num(:,1)));
        if m <= 2, kind = 'explicit'; else, kind = 'implicit'; end
        ttl = sprintf('%s (%s)', names{m}, kind);
        if mx > 3
            ex = floor(log10(mx));
            ttl = sprintf('%s, unstable: $\\max|x|\\approx %.1f\\times10^{%d}$', ...
                          ttl, mx/10^ex, ex);
        end

        plot(ax, t_num, X_num(:,1), '-o', ...
            'Color', colors(m,:), ...
            'MarkerSize', 5, ...
            'MarkerFaceColor', colors(m,:), ...
            'MarkerEdgeColor', colors(m,:), ...
            'LineWidth', 1.6, ...
            'DisplayName', label);

        xlabel(ax, 'Time, $t$ (-)', 'FontSize', 16);
        ylabel(ax, 'Solution, $x(t)$ (-)', 'FontSize', 16);
        title(ax, ttl, 'FontSize', 17, 'Interpreter', 'latex');
        legend(ax, 'Location', 'northwest', 'FontSize', 13);

        grid(ax, 'on');
        box(ax, 'on');
        xlim(ax, stability_tspan);
        ylim(ax, [-2 4.5]);          % head-room: legend sits above all data; blow-up leaves the frame
    end

    sgtitle(sprintf('Integration stability, $h_{\\mathrm{ref}}=%.2f$', stability_h), ...
            'FontSize', 20, 'Interpreter', 'latex');

    save_figure(f, fullfile(outdir, stability_names{s}));
end

%% 13. TABLES (#9, #12, #15): LaTeX files + CSV
% Compile each .tex (Overleaf, or pdflatex) to get the PDF to submit,
% or type the printed values into Word and save as PDF.

write_table(fullfile(outdir, '09_table_local.tex'), ...
    ['Estimated local truncation error scaling exponent $p$ versus step size $h$ ' ...
     '($\epsilon\approx k h^{p}$), for each method and test function'], ...
    pLocal, pDiff);

write_table(fullfile(outdir, '12_table_global_step.tex'), ...
    ['Estimated global truncation error scaling exponent $p$ versus step size $h$ ' ...
     '($\epsilon\approx k h^{p}$), for each method and test function'], ...
    pGlobal, pDiff);

write_table(fullfile(outdir, '15_table_global_calls.tex'), ...
    ['Estimated global truncation error scaling exponent $p$ versus number of rate ' ...
     'function calls $N_f$ ($\epsilon\approx k N_f^{p}$), for each method and test function'], ...
    pCalls, [NaN NaN]);

writetable(LocalErrorTable,  fullfile(outdir, 'local_error_scaling.csv'));
writetable(GlobalStepTable,  fullfile(outdir, 'global_error_step_scaling.csv'));
writetable(GlobalCallsTable, fullfile(outdir, 'global_error_call_scaling.csv'));

disp(['All deliverables saved in folder: ' outdir]);

%% LOCAL HELPER FUNCTIONS
% These helpers do not replace or modify any of your ODE methods.

function save_figure(f, base)
    try
        exportgraphics(f, [base '.png'], 'Resolution', 600);
    catch
        print(f, [base '.png'], '-dpng', '-r600');
    end
end

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

function [xf, yf, p] = fit_line(x, y, win)
    x = x(:)';
    y = y(:)';
    xf = x(win);
    coeff = polyfit(log(xf), log(y(win)), 1);
    p = coeff(1);
    yf = exp(polyval(coeff, log(xf)));
end

function write_table(fname, caption, P, pDiff)
    % Rows in the Canvas template order: FE, BE, EM, IM.
    % P rows = tests; P columns follow the script order FE, EM, BE, IM.
    row_names = {'Forward Euler', 'Backward Euler', 'Explicit Midpoint', 'Implicit Midpoint'};
    col_for_row = [1, 3, 2, 4];

    fid = fopen(fname, 'w');
    fprintf(fid, '\\documentclass[border=14pt,varwidth=17cm]{standalone}\n');
    fprintf(fid, '\\usepackage{amsmath,array}\\renewcommand{\\arraystretch}{1.4}\n');
    fprintf(fid, '\\begin{document}\n\\centering\n');
    fprintf(fid, '\\textbf{%s}\\\\[8pt]\n', caption);
    fprintf(fid, '\\begin{tabular}{|c|c|c|}\\hline\n');
    fprintf(fid, 'Method & Test 1: Estimate for $p$ (-) & Test 2: Estimate for $p$ (-) \\\\ \\hline\n');
    if any(isnan(pDiff))
        fprintf(fid, '$|X(t+h)-X(t)|$ & --- & --- \\\\ \\hline\n');
    else
        fprintf(fid, '$|X(t+h)-X(t)|$ & %.2f & %.2f \\\\ \\hline\n', pDiff(1), pDiff(2));
    end
    for r = 1:4
        c = col_for_row(r);
        fprintf(fid, '%s & %.2f & %.2f \\\\ \\hline\n', row_names{r}, P(1,c), P(2,c));
    end
    fprintf(fid, '\\end{tabular}\n\\end{document}\n');
    fclose(fid);
end
