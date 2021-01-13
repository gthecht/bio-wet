function [points, lambdas] = get_equilibriums(I)
    N = 1000;
    V_Na = 55;
    V_K  = -90;
    n0 = 0.5;

    V_vec = linspace(V_K+0.01, V_Na-0.5, N)';
    n_nullcline = nullcline_n(V_vec);

    v_nullcline = zeros(N,1);
    for ii = 1:N
        v = V_vec(ii);
        v_nullcline(ii) = fzero(@(n)V_dot_2D(v, n, I), n0);
    end
    
    null_diff = (n_nullcline - v_nullcline) .^ 2;
    [~, zero_inds] = findpeaks(-null_diff);

    V_min = V_vec(zero_inds);
    n_min = v_nullcline(zero_inds);

    delta = 1E-8;
    n_dot_dot_V = (n_dot(V_min + delta, n_min) - n_dot(V_min - delta, n_min)) / (2 * delta);
    V_dot_dot_V = (V_dot_2D(V_min + delta, n_min) - V_dot_2D(V_min - delta, n_min)) / (2 * delta);
    n_dot_dot_n = (n_dot(V_min, n_min + delta) - n_dot(V_min, n_min - delta)) / (2 * delta);
    V_dot_dot_n = (V_dot_2D(V_min, n_min + delta) - V_dot_2D(V_min, n_min - delta)) / (2 * delta);

    J3 = zeros(2,2,3);
    J3(1,1,:) = V_dot_dot_V;
    J3(1,2,:) = V_dot_dot_n;
    J3(2,1,:) = n_dot_dot_V;
    J3(2,2,:) = n_dot_dot_n;

    lambdas = zeros(2,3);
    points = zeros(2,3);
    for ii = 1:3
        lambdas(:,ii) = eig(J3(:,:,ii));
        points(:,ii) = [V_min(ii), n_min(ii)];
    end
end

