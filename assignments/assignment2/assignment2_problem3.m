% Assignment 2, Problem 3.
% Gauss-Seidel and SOR for u_xx + p(x,y)u_yy + r(x,y)u = f(x,y)
% with Dirichlet boundaries at x=0, y=0, y=1 and Neumann at x=1.

nList = [16, 32, 64];
tol = 1e-8;
maxIter = 200000;
omegaGrid = unique([1.00:0.05:1.85, 1.88:0.01:1.98]);

fprintf('\nProblem 3: Gauss-Seidel and SOR for elliptic equation\n');
fprintf('Gauss-Seidel grid refinement, tol = %.1e\n', tol);
fprintf('%8s %12s %12s %16s %12s\n', 'n', 'h', 'iterations', 'max error', 'rate');
fprintf('%8s %12s %12s %16s %12s\n', '--------', '------------', '------------', '----------------', '------------');

prevErr = NaN;
gsIterList = zeros(size(nList));
gsErrList = zeros(size(nList));
gsRateList = NaN(size(nList));
for idx = 1:numel(nList)
    n = nList(idx);
    [UGs, iterGs, x, y, errGs] = solveProblem3(n, 1.0, tol, maxIter);
    if isnan(prevErr)
        rate = NaN;
    else
        rate = log(prevErr / errGs) / log(n / nList(idx-1));
    end
    fprintf('%8d %12.5g %12d %16.8e %12.6f\n', n, 1/n, iterGs, errGs, rate);
    gsIterList(idx) = iterGs;
    gsErrList(idx) = errGs;
    gsRateList(idx) = rate;
    prevErr = errGs;
end

fprintf('\nSOR omega scan, tol = %.1e\n', tol);
fprintf('The Gauss-Seidel refinement table above is complete; the scans below only compare iteration counts.\n');

for idx = 1:numel(nList)
    n = nList(idx);
    bestIter = Inf;
    bestOmega = NaN;
    bestErr = NaN;

    fprintf('\nSOR omega scan for n = %d\n', n);
    fprintf('%10s %12s %16s %12s\n', 'omega', 'iterations', 'max error', 'converged');
    for omega = omegaGrid
        [~, iterSor, ~, ~, errSor, converged] = solveProblem3(n, omega, tol, maxIter);
        fprintf('%10.2f %12d %16.8e %12d\n', omega, iterSor, errSor, converged);
        if converged && iterSor < bestIter
            bestIter = iterSor;
            bestOmega = omega;
            bestErr = errSor;
        end
    end

    fprintf('Best SOR for n=%d: omega=%.2f, iterations=%d, max error=%.8e\n', ...
        n, bestOmega, bestIter, bestErr);

    if n == 32
        [USor, iterSor, x, y, errSor] = solveProblem3(n, bestOmega, tol, maxIter);
        [X, Y] = ndgrid(x(2:end), y(2:end-1));
        UExact = exactProblem3(X, Y);

        fig1 = figure('Name', 'Problem 3 solution n=32', 'Position', [100 100 1100 450]);
        subplot(1,2,1);
        mesh(Y, X, USor);
        title(sprintf('SOR solution, n=32, omega=%.2f', bestOmega));
        xlabel('y'); ylabel('x'); zlabel('u');
        subplot(1,2,2);
        mesh(Y, X, abs(USor - UExact));
        title(sprintf('Absolute error, iter=%d, max=%.3e', iterSor, errSor));
        xlabel('y'); ylabel('x'); zlabel('|e|');
    end
end

function [U, iter, x, y, maxErr, converged] = solveProblem3(n, omega, tol, maxIter)
    h = 1 / n;
    x = linspace(0, 1, n+1);
    y = linspace(0, 1, n+1);
    U = zeros(n, n-1);       % rows: x_1,...,x_n; cols: y_1,...,y_{n-1}
    converged = false;

    for iter = 1:maxIter
        maxChange = 0;
        for j = 1:n-1
            yj = y(j+1);
            for i = 1:n
                xi = x(i+1);
                pij = pcoef(xi, yj);
                rij = rcoef(xi, yj);

                center = -2/h^2 - 2*pij/h^2 + rij;

                if i == 1
                    leftTerm = 0;       % u(0,y)=0
                    rightTerm = U(i+1,j) / h^2;
                elseif i == n
                    leftTerm = 2 * U(i-1,j) / h^2;
                    rightTerm = 0;
                else
                    leftTerm = U(i-1,j) / h^2;
                    rightTerm = U(i+1,j) / h^2;
                end

                if j == 1
                    downTerm = 0;       % u(x,0)=0
                else
                    downTerm = pij * U(i,j-1) / h^2;
                end

                if j == n-1
                    upTerm = 0;         % u(x,1)=0
                else
                    upTerm = pij * U(i,j+1) / h^2;
                end

                rhs = fProblem3(xi, yj);
                if i == n
                    rhs = rhs - 2 * neumannProblem3(yj) / h;
                end

                gsValue = (rhs - leftTerm - rightTerm - downTerm - upTerm) / center;
                newValue = (1 - omega) * U(i,j) + omega * gsValue;
                maxChange = max(maxChange, abs(newValue - U(i,j)));
                U(i,j) = newValue;
            end
        end

        if maxChange < tol
            converged = true;
            break;
        end
    end

    [X, Y] = ndgrid(x(2:end), y(2:end-1));
    maxErr = max(max(abs(U - exactProblem3(X, Y))));
end

function val = pcoef(x, y)
    val = 1 + x.^2 + y.^2;
end

function val = rcoef(x, y)
    val = -x .* y;
end

function val = exactProblem3(x, y)
    val = sin(pi*x) .* sin(pi*y);
end

function val = neumannProblem3(y)
    val = -pi * sin(pi*y);
end

function val = fProblem3(x, y)
    u = exactProblem3(x, y);
    val = (-pi^2 * (1 + pcoef(x, y)) + rcoef(x, y)) .* u;
end
