function val = heat4_exact(t, x)
val = cos(t) .* x.^2 .* sin(pi .* x);
end

