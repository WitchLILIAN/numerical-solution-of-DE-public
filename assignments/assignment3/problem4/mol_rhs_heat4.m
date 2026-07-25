function dydt = mol_rhs_heat4(t, y, h, x)
% RHS for unknowns at x_1, ..., x_{m-1}. The Neumann boundary value at
% x=1 is eliminated using a second-order backward derivative formula.

m = numel(x) - 1;
u = reconstruct_heat4_solution(t, y, h);
dydt = zeros(m - 1, 1);

for i = 1:m - 1
    node = i + 1;
    d2u = (u(node - 1) - 2 * u(node) + u(node + 1)) / h^2;
    dydt(i) = d2u + heat4_source(t, x(node));
end

end

