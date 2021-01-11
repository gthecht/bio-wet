function n = nullcline_n(V)
% Get m = alpha(V) / (alpha(V) + beta(V))
    u = (V + 74.44);
    a = (0.1 - 0.01 * u + eps) ./ (exp(1 - 0.1 * u) - 1 + eps);
    b = 0.125 * exp(-u / 80);
    n = a ./ (a + b);
end

