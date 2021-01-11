function nv_diff = nv_dot_diff(V, n)
    ndot = n_dot(V, n);
    m = stable_m(V);
    vdot = V_dot(m, n, 1-n, V, 0);
    nv_diff = vdot - ndot;
end

