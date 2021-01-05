function [v_new, m_new, h_new, n_new] = HH_step(V, I, m, h, n, dt)
    [v_d, m_d, n_d, h_d] = HH_dot(V, I, m, h, n);
    v_new = v + v_d * dt;
    m_new = m + m_d * dt;
    h_new = h + h_d * dt;
    n_new = n + n_d * dt;
end

