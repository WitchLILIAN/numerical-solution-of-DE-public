clear; clc;

m = 40;
dt = 1 / m;
tfinal = 5.0;

[x, uNumerical, uExact, maxError, cpuTime] = backward_euler_heat5(m, dt, tfinal);

fprintf('Problem 5 backward Euler: m=%d, dt=%.5f, tfinal=%.1f\n', m, dt, tfinal);
fprintf('Max error = %.8e\n', maxError);
fprintf('CPU time  = %.6f seconds\n', cpuTime);

figure(3); clf;
plot(x, uExact, 'k-', 'LineWidth', 1.6); hold on;
plot(x, uNumerical, 'ro', 'MarkerSize', 4);
grid on;
xlabel('x'); ylabel('u(x,5)');
legend('Exact', 'Backward Euler', 'Location', 'best');
title('Problem 5 backward Euler solution, m = 40');
hold off;

figure(4); clf;
plot(x, abs(uNumerical - uExact), 'b-o', 'MarkerSize', 4);
grid on;
xlabel('x'); ylabel('absolute error at t = 5');
title('Problem 5 backward Euler error, m = 40');
