function V_dot = V_dot(m, n, h, V, I)
    %V_dot calculates V_dot
    C_m = 1;
	g_Na = 120;
	g_K = 36;
	g_CL = 0.3;

    V_Na = 55;
	V_K  = -90;
	V_CL = -60;

	vna = g_Na * (m .^ 3) .* h .* (V - V_Na);
	vk  = g_K * (n .^ 4) .* (V - V_K);
	vcl = g_CL * (V - V_CL);
	V_dot = - (1 / C_m) * (vna + vk + vcl) + I / C_m;
end