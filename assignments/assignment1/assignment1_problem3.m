clear; close all; clc

% Adjustable parameters in the assignment statement.
K = 1;
theta1 = pi/6;
theta2 = pi/3;
tol = 1.e-10;
kmax = 50;

mList = [40 80 160]';

fprintf('\nProblem 3: nonlinear pendulum BVP\n')
fprintf('K = %.8g, theta1 = %.8g, theta2 = %.8g\n', K, theta1, theta2)
fprintf('%8s %14s %12s %18s %18s %12s\n', ...
    'm', 'h', 'Iter', 'ResidualInf', 'CorrectionInf', 'Converged')
for j = 1:length(mList)
    [t,theta,info] = pendulum_newton(K,theta1,theta2,mList(j),tol,kmax);
    fprintf('%8d %14.8e %12d %18.8e %18.8e %12d\n', ...
        mList(j), info.h, info.iterations, info.residualNorm, ...
        info.correctionNorm, info.converged)
end

[t,theta,info] = pendulum_newton(K,theta1,theta2,80,tol,kmax);

figure(1)
plot(t, theta, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4)
grid on
xlabel('t')
ylabel('\theta(t)')
title(sprintf('K=%g, theta_1=%.4f, theta_2=%.4f', K, theta1, theta2))

fprintf('\nSample m=80 solve info:\n')
fprintf('iterations     = %d\n', info.iterations)
fprintf('correction inf = %.8e\n', info.correctionNorm)
fprintf('residual inf   = %.8e\n', info.residualNorm)
