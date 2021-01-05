function [m_d, n_d, h_d, v_d] = HH_dot(m, n, h, V, I)
    m_d = m_dot(V, m);
    n_d = n_dot(V, n);
    h_d = h_dot(V, h);
    v_d = V_dot(m, n, h, V, I);
end

