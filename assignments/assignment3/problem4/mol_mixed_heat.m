function [x, uFull, uExact, maxError, cpuTime] = mol_mixed_heat(m, tfinal)
% Method of lines with ode15s for the same spatial discretization family.

h = 1 / m;
x = linspace(0, 1, m + 1);
y0 = heat4_exact(0, x(2:m)).';

opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
timer = tic;
[~, y] = ode15s(@(t, y) mol_rhs_heat4(t, y, h, x), [0, tfinal], y0, opts);
cpuTime = toc(timer);

yFinal = y(end, :).';
uFull = reconstruct_heat4_solution(tfinal, yFinal, h);
uExact = heat4_exact(tfinal, x);
maxError = max(abs(uFull - uExact));

end

