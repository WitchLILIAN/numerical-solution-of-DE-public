function u = reconstruct_heat4_solution(t, y, h)
% Build full grid values from interior unknowns and mixed boundary data.

m = numel(y) + 1;
u = zeros(1, m + 1);
u(1) = heat4_g1(t);
u(2:m) = y(:).';
u(m + 1) = (4 * u(m) - u(m - 1) + 2 * h * heat4_g2(t)) / 3;

end

