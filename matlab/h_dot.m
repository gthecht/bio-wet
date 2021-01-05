function h_dot = h_dot(V, h)
    u = (V + 74.44);
    a = 0.07 * exp(-u / 20);
    b = 1 / (exp(3 - 0.1 * u) + 1 + eps);
    h_dot = get_dot(a, b, h);
end
