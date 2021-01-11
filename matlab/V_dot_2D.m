function v_dot = V_dot_2D(V, n, I)
    if nargin < 3
        I = 0;
    end
    m = stable_m(V);
    v_dot = V_dot(m, n, 1-n, V, I);
end

