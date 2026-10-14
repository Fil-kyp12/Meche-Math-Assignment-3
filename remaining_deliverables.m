
%% REMAINING DELIVERABLES
% Uses existing integration methods, rate functions, and exact solutions.
% Does not modify any existing functions.

clearvars -except forward_euler_step explicit_midpoint_step ...
    backward_euler_step implicit_midpoint_step ...
    fixed_step_integration rate_func01 rate_func02 ...
    solution01 solution02

close all;
clc;

%% Method definitions

methods = {@forward_euler_step, @explicit_midpoint_step, ...
           @backward_euler_step, @implicit_midpoint_step};

names = {'Forward Euler', 'Explicit Midpoint', ...
         'Backward Euler', 'Implicit Midpoint'};

nMethods = length(methods);

%% Test problem definitions

rates = {@rate_func01, @rate_func02};
solutions = {@solution01, @solution02};

tspan = [0, 10];
tref = 0.5;       % Nonzero reference time for local error

X0 = {solution01(0), solution02(0)};

% Candidate step sizes, spanning multiple orders of magnitude.
h_values = logspace(log10(0.005), log10(0.5), 10);
nH = length(h_values);

%% Storage arrays
% Dimensions: method x step size x test problem

localErr  = zeros(nMethods, nH, 2);
globalErr = zeros(nMethods, nH, 2);
calls     = zeros(nMethods, nH, 2);
hActual   = zeros(nMethods, nH, 2);

%% Calculate local and global errors for both test problems

for test = 1:2

    rate = rates{test};
    sol  = solutions{test};

    % Local errors: one numerical step starting from the exact solution.
    Xref = sol(tref);

    for m = 1:nMethods
        for k = 1:nH

            h = h_values(k);

            [Xnext, ~] = methods{m}(rate, tref, Xref, h);

            XexactNext = sol(tref + h);

            localErr(m,k,test) = norm( ...
                Xnext(:) - XexactNext(:), 2);

        end
    end

    % Global errors: integrate over the same IVP and interval.
    for m = 1:nMethods
        for k = 1:nH

            [t, X, havg, neval] = fixed_step_integration( ...
                rate, methods{m}, tspan, X0{test}, h_values(k));

            XexactEnd = sol(t(end));
            XnumericalEnd = X(end,:)';

            globalErr(m,k,test) = norm( ...
                XnumericalEnd(:) - XexactEnd(:), 2);

            calls(m,k,test) = neval;
            hActual(m,k,test) = havg;

        end
    end
end

%% TABLE 1: Local truncation error scaling exponent p

pLocal = zeros(nMethods,2);

for test = 1:2
    for m = 1:nMethods
        pLocal(m,test) = estimate_slope( ...
            h_values, squeeze(localErr(m,:,test)));
    end
end

LocalErrorScalingTable = table( ...
    names', pLocal(:,1), pLocal(:,2), ...
    'VariableNames', {'Method','Test1_p','Test2_p'});

disp('LOCAL TRUNCATION ERROR SCALING');
disp(LocalErrorScalingTable);

%% TABLE 2: Global truncation error scaling with step size

pGlobalH = zeros(nMethods,2);

for test = 1:2
    for m = 1:nMethods
        pGlobalH(m,test) = estimate_slope( ...
            squeeze(hActual(m,:,test)), ...
            squeeze(globalErr(m,:,test)));
    end
end

GlobalErrorStepTable = table( ...
    names', pGlobalH(:,1), pGlobalH(:,2), ...
    'VariableNames', {'Method','Test1_p','Test2_p'});

disp('GLOBAL ERROR SCALING WITH STEP SIZE');
disp(GlobalErrorStepTable);

%% TABLE 3: Global truncation error scaling with rate-function calls

pCalls = zeros(nMethods,2);

for test = 1:2
    for m = 1:nMethods
        pCalls(m,test) = estimate_slope( ...
            squeeze(calls(m,:,test)), ...
            squeeze(globalErr(m,:,test)));
    end
end

GlobalErrorCallsTable = table( ...
    names', pCalls(:,1), pCalls(:,2), ...
    'VariableNames', {'Method','Test1_p','Test2_p'});

disp('GLOBAL ERROR SCALING WITH RATE-FUNCTION CALLS');
disp(GlobalErrorCallsTable);

%% PLOT 1: Local truncation error vs. step size
% All four methods. No fit lines or analytical-difference curve.
% Uses Test 1 and nonzero tref.

figure;
hold on;

for m = 1:nMethods
    loglog(h_values, squeeze(localErr(m,:,1)), 'o-', ...
        'LineWidth', 1.5, 'DisplayName', names{m});
end

xlabel('Step size h');
ylabel('Local truncation error');
title(sprintf('Local Error vs. Step Size (t_{ref} = %.2f)', tref));
legend('Location','best');
grid on;
hold off;

%% PLOT 2: Local error, explicit methods, with fits
% Includes |X(tref+h)-X(tref)| for comparison.

figure;
hold on;

explicitIdx = [1, 2];

for j = 1:length(explicitIdx)
    m = explicitIdx(j);
    E = squeeze(localErr(m,:,1));

    loglog(h_values, E, 'o', 'LineWidth', 1.5, ...
        'DisplayName', [names{m} ' data']);

    p = estimate_slope(h_values, E);
    fitLine = fitted_values(h_values, E);

    loglog(h_values, fitLine, '--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s fit, p = %.2f', names{m}, p));
end

% Analytical difference, as requested for this plot only.
analyticalDifference = zeros(size(h_values));

for k = 1:nH
    analyticalDifference(k) = abs( ...
        solution01(tref + h_values(k)) - solution01(tref));
end

loglog(h_values, analyticalDifference, 'k-.', ...
    'LineWidth', 1.5, ...
    'DisplayName', '|X(tref+h)-X(tref)|');

xlabel('Step size h');
ylabel('Error / analytical difference');
title('Local Error Scaling: Explicit Methods');
legend('Location','best');
grid on;
hold off;

%% PLOT 3: Global truncation error vs. step size
% All four methods. No fit lines.

figure;
hold on;

for m = 1:nMethods
    loglog(squeeze(hActual(m,:,1)), ...
        squeeze(globalErr(m,:,1)), 'o-', ...
        'LineWidth', 1.5, 'DisplayName', names{m});
end

xlabel('Step size h');
ylabel('Global endpoint error');
title('Global Error vs. Step Size');
legend('Location','best');
grid on;
hold off;

%% PLOT 4: Global error vs. step size, explicit methods with fits

figure;
hold on;

for j = 1:length(explicitIdx)
    m = explicitIdx(j);

    h = squeeze(hActual(m,:,1));
    E = squeeze(globalErr(m,:,1));

    loglog(h, E, 'o', 'LineWidth', 1.5, ...
        'DisplayName', [names{m} ' data']);

    p = estimate_slope(h, E);
    loglog(h, fitted_values(h, E), '--', ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s fit, p = %.2f', names{m}, p));
end

xlabel('Step size h');
ylabel('Global endpoint error');
title('Global Error vs. Step Size: Explicit Methods');
legend('Location','best');
grid on;
hold off;

%% PLOT 5: Global error vs. rate-function calls
% All four methods. No fit lines.

figure;
hold on;

for m = 1:nMethods
    loglog(squeeze(calls(m,:,1)), ...
        squeeze(globalErr(m,:,1)), 'o-', ...
        'LineWidth', 1.5, 'DisplayName', names{m});
end

xlabel('Number of rate-function calls');
ylabel('Global endpoint error');
title('Global Error vs. Rate-Function Calls');
legend('Location','best');
grid on;
hold off;

%% PLOT 6: Global error vs. rate-function calls
% Explicit methods with fitted lines.

figure;
hold on;

for j = 1:length(explicitIdx)
    m = explicitIdx(j);

    Nf = squeeze(calls(m,:,1));
    E = squeeze(globalErr(m,:,1));

    loglog(Nf, E, 'o', 'LineWidth', 1.5, ...
        'DisplayName', [names{m} ' data']);

    p = estimate_slope(Nf, E);
    loglog(Nf, fitted_values(Nf, E), '--', ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s fit, p = %.2f', names{m}, p));
end

xlabel('Number of rate-function calls');
ylabel('Global endpoint error');
title('Global Error vs. Rate-Function Calls: Explicit Methods');
legend('Location','best');
grid on;
hold off;

%% PLOT 7: Stability experiment, href = 0.38
% Test 1, interval [0,20], four stacked subplots.

make_stability_figure(0.38, methods, names);

%% PLOT 8: Stability experiment, href = 0.45

make_stability_figure(0.45, methods, names);

%% Optional: save the three tables as CSV files.

writetable(LocalErrorScalingTable, 'local_error_scaling.csv');
writetable(GlobalErrorStepTable, 'global_error_step_scaling.csv');
writetable(GlobalErrorCallsTable, 'global_error_call_scaling.csv');

disp('All requested remaining tables and plots have been generated.');

%% LOCAL HELPER FUNCTIONS
% These helpers are added only in this script.
% Your existing functions remain unchanged.

function p = estimate_slope(x, y)

    x = x(:);
    y = y(:);

    valid = isfinite(x) & isfinite(y) & x > 0 & y > 0;

    if sum(valid) < 2
        p = NaN;
        return;
    end

    coefficients = polyfit(log(x(valid)), log(y(valid)), 1);
    p = coefficients(1);

end

function yfit = fitted_values(x, y)

    x = x(:);
    y = y(:);

    valid = isfinite(x) & isfinite(y) & x > 0 & y > 0;

    coefficients = polyfit(log(x(valid)), log(y(valid)), 1);

    yfit = exp(polyval(coefficients, log(x)));

end

function make_stability_figure(href, methods, names)

    tspan = [0, 20];
    X0 = solution01(0);

    tFine = linspace(tspan(1), tspan(2), 3000);
    XFine = solution01(tFine);

    figure('Name', sprintf('Stability href = %.2f', href));

    for m = 1:length(methods)

        [t, X, havg] = fixed_step_integration( ...
            @rate_func01, methods{m}, tspan, X0, href);

        subplot(4,1,m);

        plot(tFine, XFine, 'k-', 'LineWidth', 1.5, ...
            'DisplayName', 'Analytical solution');
        hold on;

        plot(t, X(:,1), 'o-', 'MarkerSize', 3, ...
            'DisplayName', 'Numerical solution');

        ylabel('x(t)');
        title(sprintf('%s, h_{avg} = %.4f', names{m}, havg));
        legend('Location','best');
        grid on;
        hold off;

    end

    xlabel('Time t');

    sgtitle(sprintf('Integration Stability: h_{ref} = %.2f', href));

end
