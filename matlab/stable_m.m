function m = stable_m(V)
% Get m = alpha(V) / (alpha(V) + beta(V))
    u = (V + 74.44);
    a = (2.5 - 0.1 * u + eps) ./ (exp(2.5 - 0.1 * u) - 1 + eps);
    b = 4 .* exp(-u / 18);
    m = a ./ (a + b);
end

