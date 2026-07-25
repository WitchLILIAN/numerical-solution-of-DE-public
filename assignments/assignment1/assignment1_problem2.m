clear; close all; clc

nList = [32 64 128 256 512 1024]';
K = 1;
a = 0;
b = 1;

% Problem 2(a): u(x) = sin(5x), so f(x) = u''(x) - u(x) = -26 sin(5x).
exactA = @(x) sin(5*x);
fA = @(x) -26*sin(5*x);
alphaA = exactA(a);
betaA = exactA(b);
[hA, errA, rateA, lastXA, lastUA] = run_case(a,b,alphaA,betaA,fA,exactA,nList,K);
fprintf('\nProblem 2(a): u(x)=sin(5x)\n')
fprintf('f(x) = -26 sin(5x), alpha = %.16g, beta = %.16g\n', alphaA, betaA)
print_error_table(nList,hA,errA,rateA)

figure(1)
loglog(hA, errA, 'o-', 'LineWidth', 1.5, 'MarkerSize', 7)
grid on
xlabel('h')
ylabel('maximum error')
title('Problem 2(a): u(x)=sin(5x)')

figure(2)
plot(lastXA, lastUA, 'o', lastXA, exactA(lastXA), '-', ...
    'LineWidth', 1.4, 'MarkerSize', 4)
grid on
xlabel('x')
ylabel('u(x)')
legend('finite difference', 'exact', 'Location', 'best')
title('Problem 2(a): numerical and exact solution, n=1024')

% Problem 2(b): f(x)=x^3 and u(0)=u(1)=0.
% Analytic solution: u(x)=7 sinh(x)/sinh(1) - x^3 - 6x.
exactB = @(x) 7*exp(1)*exp(x)./(exp(2)-1)+7*exp(1)*exp(-x)./(1-exp(2)) - x.^3 - 6*x;
fB = @(x) x.^3;
alphaB = 0;
betaB = 0;
[hB, errB, rateB, lastXB, lastUB] = run_case(a,b,alphaB,betaB,fB,exactB,nList,K);
fprintf('\nProblem 2(b): f(x)=x^3\n')
fprintf('exact solution: u(x) = 7*sinh(x)/sinh(1) - x^3 - 6*x\n')
fprintf('alpha = %.16g, beta = %.16g\n', alphaB, betaB)
print_error_table(nList,hB,errB,rateB)

figure(3)
loglog(hB, errB, 'o-', 'LineWidth', 1.5, 'MarkerSize', 7)
grid on
xlabel('h')
ylabel('maximum error')
title('Problem 2(b): f(x)=x^3')

figure(4)
plot(lastXB, lastUB, 'o', lastXB, exactB(lastXB), '-', ...
    'LineWidth', 1.4, 'MarkerSize', 4)
grid on
xlabel('x')
ylabel('u(x)')
legend('finite difference', 'exact', 'Location', 'best')
title('Problem 2(b): numerical and exact solution, n=1024')

function [h, err, rate, lastX, lastU] = run_case(a,b,ua,ub,f,exact,nList,K)
    err = zeros(size(nList));
    h = zeros(size(nList));
    for j = 1:length(nList)
        n = nList(j);
        h(j) = (b-a)/n;
        [x,U] = two_point_reaction(a,b,ua,ub,f,n,K);
        err(j) = norm(U - exact(x), inf);
        lastX = x;
        lastU = U;
    end

    rate = nan(size(err));
    for j = 2:length(err)
        rate(j) = log(err(j-1)/err(j))/log(2);
    end
end

function print_error_table(nList,h,err,rate)
    fprintf('%8s %14s %18s %12s\n', 'n', 'h', 'MaxError', 'Rate')
    for j = 1:length(nList)
        if j == 1
            fprintf('%8d %14.8e %18.8e %12s\n', nList(j), h(j), err(j), '-')
        else
            fprintf('%8d %14.8e %18.8e %12.6f\n', nList(j), h(j), err(j), rate(j))
        end
    end
end
