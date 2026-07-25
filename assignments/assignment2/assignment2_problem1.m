% Assignment 2, Problem 1.
% Compare central and central-upwind finite differences for
% eps*u'' - u' = -1, 0<x<1, u(0)=1, u(1)=3.
clear; clc;

epsList = [0.3, 0.1, 0.05, 0.0005];
nRefDefault = [10, 25, 50, 100, 200, 400, 800];
nRefLayer = [10, 25, 50, 100, 200, 400, 800, 1600, 3200, 6400];
nPlot = [10, 25, 100];
methodNames = {'central', 'upwind'};

fprintf('\nProblem 1: grid refinement for convection-diffusion equation\n');
fprintf('Exact solution: u(x)=1+x+(exp(x/eps)-1)/(exp(1/eps)-1)\n');
for ie = 1:numel(epsList)
    epsVal = epsList(ie);
    if epsVal <= 1e-3
        nRef = nRefLayer;
    else
        nRef = nRefDefault;
    end

    fig = figure('Name', sprintf('Problem 1 epsilon=%g', epsVal), ...
        'Position', [100 100 1150 700]);

    for im = 1:numel(methodNames)
        method = methodNames{im};
        prevErr = NaN;

        fprintf('\nepsilon = %.4g, method = %s\n', epsVal, method);
        fprintf('%8s %12s %16s %12s\n', 'n', 'h', 'max error', 'rate');
        

        for in = 1:numel(nRef)
            n = nRef(in);
            [x, uNum] = solveConvectionDiffusion(epsVal, n, method);
            uEx = exactConvectionDiffusion(x, epsVal);
            err = max(abs(uNum - uEx));
            if isnan(prevErr)
                rate = NaN;
            else
                rate = log(prevErr / err) / log(n / nRef(in-1));
            end
            fprintf('%8d %12.5g %16.8e %12.6f\n', n, 1/n, err, rate);
            prevErr = err;
        end

        for ip = 1:numel(nPlot)
            n = nPlot(ip);
            [x, uNum] = solveConvectionDiffusion(epsVal, n, method);
            uEx = exactConvectionDiffusion(x, epsVal);
            subplot(numel(methodNames), numel(nPlot), (im-1)*numel(nPlot)+ip);
            plot(x, uEx, 'k-', 'LineWidth', 1.4); hold on;
            plot(x, uNum, 'ro-', 'MarkerSize', 3, 'LineWidth', 1.0);
            grid on;
            title(sprintf('%s, h=%.4g', method, 1/n), 'Interpreter', 'none');
            xlabel('x'); ylabel('u');
            legend('exact', 'computed', 'Location', 'best');
        end
    end

    sgtitle(sprintf('Problem 1: epsilon = %.4g', epsVal));
end

function [x, U] = solveConvectionDiffusion(epsVal, n, method)
    h = 1 / n;
    x = linspace(0, 1, n+1)';
    N = n - 1;
    rhs = -ones(N, 1);

    switch method
        case 'central'
            lower = epsVal / h^2 + 1 / (2*h);
            diagv = -2 * epsVal / h^2;
            upper = epsVal / h^2 - 1 / (2*h);
        case 'upwind'
            lower = epsVal / h^2 + 1 / h;
            diagv = -2 * epsVal / h^2 - 1 / h;
            upper = epsVal / h^2;
        otherwise
            error('Unknown method: %s', method);
    end

    A = spdiags([lower*ones(N,1), diagv*ones(N,1), upper*ones(N,1)], ...
        -1:1, N, N);
    rhs(1) = rhs(1) - lower * 1;
    rhs(end) = rhs(end) - upper * 3;

    U = zeros(n+1, 1);
    U(1) = 1;
    U(end) = 3;
    U(2:n) = A \ rhs;
end

function u = exactConvectionDiffusion(x, epsVal)
    % Stable form of (exp(x/eps)-1)/(exp(1/eps)-1), avoiding overflow.
    den = 1 - exp(-1/epsVal);
    layer = (exp((x - 1) / epsVal) - exp(-1/epsVal)) / den;
    u = 1 + x + layer;
end
