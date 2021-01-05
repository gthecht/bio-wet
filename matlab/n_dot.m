function n_dot = n_dot(V, n)
    u = (V + 74.44);
    a = (0.1 - 0.01 * u + eps) / (exp(1 - 0.1 * u) - 1 + eps);
    b = 0.125 * exp(-u / 80);
    n_dot = get_dot(a, b, n);
end
