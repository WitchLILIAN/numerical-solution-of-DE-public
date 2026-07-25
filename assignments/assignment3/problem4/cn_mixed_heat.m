function [x, uFull, uExact, maxError, cpuTime] = cn_mixed_heat(m, dt, tfinal)
% Crank-Nicolson method for u_t = u_xx + f on 0 < x < 1.
% Boundary conditions: u(0,t)=g1(t), u_x(1,t)=g2(t).

h = 1 / m;
x = linspace(0, 1, m + 1);
r = dt / (2 * h^2);
numUnknowns = m;

A = sparse(numUnknowns, numUnknowns);
B = sparse(numUnknowns, numUnknowns);

for i = 1:m - 1
    A(i, i) = 1 + 2 * r;
    B(i, i) = 1 - 2 * r;
    if i > 1
        A(i, i - 1) = -r;
        B(i, i - 1) = r;
    end
    if i < m
        A(i, i + 1) = -r;
        B(i, i + 1) = r;
    end
end

% Second-order Neumann boundary at x=1:
% (3 U_m - 4 U_{m-1} + U_{m-2}) / (2h) = g2(t).
A(m, m) = 1;
A(m, m - 1) = -4 / 3;
A(m, m - 2) = 1 / 3;

u = heat4_exact(0, x(2:end)).';
numSteps = round(tfinal / dt);
t = 0;

timer = tic;
for n = 1:numSteps
    tNext = t + dt;
    tMid = t + 0.5 * dt;
    rhs = B * u;

    rhs(1) = rhs(1) + r * (heat4_g1(t) + heat4_g1(tNext));
    rhs(1:m - 1) = rhs(1:m - 1) + dt * heat4_source(tMid, x(2:m)).';
    rhs(m) = (2 * h / 3) * heat4_g2(tNext);

    u = A \ rhs;
    t = tNext;
end
cpuTime = toc(timer);

uFull = [heat4_g1(t), u.'];
uExact = heat4_exact(tfinal, x);
maxError = max(abs(uFull - uExact));

end

