function dy2dt = hhx(t, y, a)
    I = a;
    [m_d, n_d, h_d, v_d] = HH_dot(y(1), y(2), y(3), y(4), I);
    dy2dt = [m_d; n_d; h_d; v_d];
end

