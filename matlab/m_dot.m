function m_dot = m_dot(V, m)
    u = (V + 74.44);
    a = (2.5 - 0.1 * u + eps) / (exp(2.5 - 0.1 * u) - 1 + eps);
    b = 4 * exp(-u / 18);
    m_dot = get_dot(a, b, m);
end
