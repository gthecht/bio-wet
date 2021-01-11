function dy2dt = hh2d(t, y, I)
    V = y(2);
    m = stable_m(V);
    
    [m_d, n_d, h_d, v_d] = HH_dot(m, y(1), 1 - y(1), y(2), I);
    dy2dt = [n_d; v_d];
end

