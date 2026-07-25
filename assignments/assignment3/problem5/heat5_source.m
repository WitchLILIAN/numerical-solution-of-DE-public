function val = heat5_source(t, x)
val = -sin(t) .* sin(pi .* x) + cos(t) .* pi^2 .* sin(pi .* x);
end

