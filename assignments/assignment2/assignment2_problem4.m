% Assignment 2, Problem 4.
% Direct sparse solve for div(p grad u) - q u = f with Dirichlet data
% taken from the analytic solution on all four boundaries.

cases = {'poly', 'trigexp'};
nList = [16, 32, 64, 128];

fprintf('\nProblem 4: self-adjoint elliptic equation\n');
for ic = 1:numel(cases)
    caseName = cases{ic};
    prevErr = NaN;
    fprintf('\ncase = %s\n', caseName);
    fprintf('%8s %12s %16s %18s %12s\n', 'n', 'h', 'max abs error', 'max relative error', 'rate');
    for in = 1:numel(nList)
        n = nList(in);
        [U, UExact, x, y, absErr, relErr] = solveSelfAdjoint(caseName, n);
        if isnan(prevErr)
            rate = NaN;
        else
            rate = log(prevErr / absErr) / log(n / nList(in-1));
        end
        fprintf('%8d %12.5g %16.8e %18.8e %12.6f\n', n, 1/n, absErr, relErr, rate);
        prevErr = absErr;

        if n == 32
            [X, Y] = ndgrid(x(2:end-1), y(2:end-1));
            pointRel = abs(U - UExact) ./ max(abs(UExact), 1e-14);
            fig = figure('Name', sprintf('Problem 4 %s n=32', caseName), ...
                'Position', [100 100 1200 420]);
            subplot(1,3,1);
            mesh(Y, X, U);
            title(sprintf('%s computed solution', caseName), 'Interpreter', 'none');
            xlabel('y'); ylabel('x'); zlabel('u');
            subplot(1,3,2);
            mesh(Y, X, abs(U - UExact));
            title('Absolute error');
            xlabel('y'); ylabel('x'); zlabel('|e|');
            subplot(1,3,3);
            mesh(Y, X, pointRel);
            title('Pointwise relative error');
            xlabel('y'); ylabel('x'); zlabel('relative error');
        end
    end
end

function [U, UExact, x, y, maxAbsErr, maxRelErr] = solveSelfAdjoint(caseName, n)
    a = 0; b = 1; c = 0; d = 1;
    hx = (b-a) / n;
    hy = (d-c) / n;
    x = linspace(a, b, n+1);
    y = linspace(c, d, n+1);

    N = (n-1) * (n-1);
    A = sparse(N, N);
    rhs = zeros(N, 1);

    for j = 1:n-1
        yj = y(j+1);
        for i = 1:n-1
            xi = x(i+1);
            k = idx(i, j, n);

            pe = pSelf(caseName, xi + hx/2, yj);
            pw = pSelf(caseName, xi - hx/2, yj);
            pn = pSelf(caseName, xi, yj + hy/2);
            ps = pSelf(caseName, xi, yj - hy/2);
            qv = qSelf(caseName, xi, yj);

            ce = pe / hx^2;
            cw = pw / hx^2;
            cn = pn / hy^2;
            cs = ps / hy^2;
            cc = -(pe + pw) / hx^2 - (pn + ps) / hy^2 - qv;

            A(k,k) = cc;
            rhs(k) = fSelf(caseName, xi, yj);

            if i == 1
                rhs(k) = rhs(k) - cw * uSelf(caseName, x(1), yj);
            else
                A(k, idx(i-1, j, n)) = cw;
            end

            if i == n-1
                rhs(k) = rhs(k) - ce * uSelf(caseName, x(end), yj);
            else
                A(k, idx(i+1, j, n)) = ce;
            end

            if j == 1
                rhs(k) = rhs(k) - cs * uSelf(caseName, xi, y(1));
            else
                A(k, idx(i, j-1, n)) = cs;
            end

            if j == n-1
                rhs(k) = rhs(k) - cn * uSelf(caseName, xi, y(end));
            else
                A(k, idx(i, j+1, n)) = cn;
            end
        end
    end

    vecU = A \ rhs;
    U = reshape(vecU, n-1, n-1);
    [X, Y] = ndgrid(x(2:end-1), y(2:end-1));
    UExact = uSelf(caseName, X, Y);
    maxAbsErr = max(max(abs(U - UExact)));
    maxRelErr = maxAbsErr / max(max(abs(UExact)));
end

function k = idx(i, j, n)
    k = i + (j-1) * (n-1);
end

function val = uSelf(caseName, x, y)
    switch caseName
        case 'poly'
            val = x.^2 + y.^2;
        case 'trigexp'
            val = cos(x) .* sin(y);
        otherwise
            error('Unknown case: %s', caseName);
    end
end

function val = pSelf(caseName, x, y)
    switch caseName
        case 'poly'
            val = ones(size(x + y));
        case 'trigexp'
            val = exp(x + y);
        otherwise
            error('Unknown case: %s', caseName);
    end
end

function val = qSelf(caseName, x, y)
    switch caseName
        case 'poly'
            val = ones(size(x + y));
        case 'trigexp'
            val = x.^2 + y.^2;
        otherwise
            error('Unknown case: %s', caseName);
    end
end

function val = fSelf(caseName, x, y)
    switch caseName
        case 'poly'
            val = 4 - x.^2 - y.^2;
        case 'trigexp'
            p = exp(x + y);
            u = cos(x) .* sin(y);
            divPart = p .* (-sin(x).*sin(y) - 2*cos(x).*sin(y) + cos(x).*cos(y));
            val = divPart - (x.^2 + y.^2) .* u;
        otherwise
            error('Unknown case: %s', caseName);
    end
end
