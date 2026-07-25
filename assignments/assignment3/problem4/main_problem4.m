clear; close all; clc;

tfinal = 5.0;
gridList = [20, 40, 80, 160];
numCases = numel(gridList);

cnErr = zeros(numCases, 1);
molErr = zeros(numCases, 1);
cnCpu = zeros(numCases, 1);
molCpu = zeros(numCases, 1);
cnRate = NaN(numCases, 1);
molRate = NaN(numCases, 1);

for q = 1:numCases
    m = gridList(q);
    dt = 1 / m;

    [xCn, uCn, uExactCn, cnErr(q), cnCpu(q)] = cn_mixed_heat(m, dt, tfinal);
    [xMol, uMol, uExactMol, molErr(q), molCpu(q)] = mol_mixed_heat(m, tfinal);

    if q > 1
        cnRate(q) = log(cnErr(q - 1) / cnErr(q)) / log(2);
        molRate(q) = log(molErr(q - 1) / molErr(q)) / log(2);
    end

    if m == 80
        xPlot = xCn;
        uCnPlot = uCn;
        uMolPlot = interp1(xMol, uMol, xPlot, 'linear');
        uExactPlot = uExactCn;
    end
end

fprintf('\nProblem 4: Crank-Nicolson and MOL for u_t = u_xx + f(x,t)\n');
fprintf('Exact solution: u(x,t) = cos(t) x^2 sin(pi x), tfinal = %.1f\n\n', tfinal);

fprintf('\n%6s %10s %16s %10s %12s %16s %10s %12s\n', ...
    'm', 'dt_CN', 'CN_MaxError', 'CN_Rate', 'CN_CPU', ...
    'MOL_MaxError', 'MOL_Rate', 'MOL_CPU');
for q = 1:numCases
    fprintf('%6d %10.5f %16.8e %10.4f %12.6f %16.8e %10.4f %12.6f\n', ...
        gridList(q), 1 / gridList(q), cnErr(q), cnRate(q), cnCpu(q), ...
        molErr(q), molRate(q), molCpu(q));
end

figure(1); clf;
plot(xPlot, uExactPlot, 'k-', 'LineWidth', 1.6); hold on;
plot(xPlot, uCnPlot, 'ro', 'MarkerSize', 4);
plot(xPlot, uMolPlot, 'b+', 'MarkerSize', 4);
grid on;
xlabel('x'); ylabel('u(x,5)');
legend('Exact', 'Crank-Nicolson', 'MOL ode15s', 'Location', 'best');
title('Problem 4 solution comparison, m = 80');
hold off;

figure(2); clf;
plot(xPlot, abs(uCnPlot - uExactPlot), 'r-o', 'MarkerSize', 4); hold on;
plot(xPlot, abs(uMolPlot - uExactPlot), 'b-+', 'MarkerSize', 4);
grid on;
xlabel('x'); ylabel('absolute error at t = 5');
legend('Crank-Nicolson', 'MOL ode15s', 'Location', 'best');
title('Problem 4 error comparison, m = 80');
hold off;
