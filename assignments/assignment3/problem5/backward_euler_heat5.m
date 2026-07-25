function [x, uFull, uExact, maxError, cpuTime] = backward_euler_heat5(m, dt, tfinal)
% Backward Euler method for u_t = u_xx + f with homogeneous Dirichlet BC.

h = 1 / m;
x = linspace(0, 1, m + 1);
interiorX = x(2:m);
numInterior = m - 1;
r = dt / h^2;

mainDiag = (1 + 2 * r) * ones(numInterior, 1);
offDiag = -r * ones(numInterior, 1);
A = spdiags([offDiag, mainDiag, offDiag], -1:1, numInterior, numInterior);

u = heat5_exact(0, interiorX).';
numSteps = round(tfinal / dt);
t = 0;

timer = tic;
for n = 1:numSteps
    t = t + dt;
    rhs = u + dt * heat5_source(t, interiorX).';
    u = A \ rhs;
end
cpuTime = toc(timer);

uFull = [0, u.', 0];
uExact = heat5_exact(tfinal, x);
maxError = max(abs(uFull - uExact));

end

