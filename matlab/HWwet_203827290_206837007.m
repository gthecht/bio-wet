close all

%% Constants:
x = [2 0 6 8 3 7 0 0 7];
y = [2 0 3 8 2 7 2 9 0];
%% PART I

%% 1.1

%% 1.2

%% 1.3

%% 2.1

%% 2.2

%% 3.1

%% 3.2

%% PART II

%% 1.1
tspan = [0, 200];
y0 = [0.05 ,0.32, 0.6, -74];

a1 = (1 + x(9) / 9);

[t,y] = ode15s(@(t,y)hhx(t,y,a1), tspan, y0);

show_part2_Q1(t, y, a1)

%{
We can see no PPs are created by this current. This shows that the system
got accostumed to the input.
%}
%% 1.2
a2 = (5 + x(8));

[t,y] = ode15s(@(t,y)hhx(t,y,a2), tspan, y0);

show_part2_Q1(t, y, a2)
%% 1.3
%{
We can see in section 1.2 that PPs are created all the time.
%}
N = 100;
a_vec = linspace(a1, a2, N);
for a = a_vec
    [t,y] = ode15s(@(t,y)hhx(t,y,a), tspan, y0);
    v = y(:,4);
    v80 = prctile(v, 80);
    last_max = find(v > v80, 1, "last") / length(v);
    if (0.6 < last_max) && (last_max < 0.7)
        a3 = a;
        break
    end
end

show_part2_Q1(t, y, a3)

% Find minimal alpha and maximal alpha:

%% 1.4
N = 1E3;
a_vec = linspace(a1, a2, N);
a_min = 0;
a_max = 0;
for a = a_vec
    [t,y] = ode15s(@(t,y)hhx(t,y,a), tspan, y0);
    v = y(:,4);
    [peaks, locs] = findpeaks(v,t);
    if (peaks(1) < 0)
        a_max = a;
    elseif (peaks(end) > 0 && a_min == 0)
        a_min = a;
    end
end
%%
figure(5);
[t,y] = ode15s(@(t,y)hhx(t,y,a_max), tspan, y0);
v = y(:,4);
p1 = subplot(1,2,1);
plot(t, v);
title("$V(t)$ for $I = " + num2str(a_max) + " \cdot u(t)$", ...
    "interpreter", "latex");
xlabel("t[msec]");
ylabel("V[mV]");

[t,y] = ode15s(@(t,y)hhx(t,y,a_min), tspan, y0);
v = y(:,4);
p2 = subplot(1,2,2);
plot(t, v);
title("$V(t)$ for $I = " + num2str(a_min) + " \cdot u(t)$", ...
    "interpreter", "latex");
xlabel("t[msec]");
ylabel("V[mV]");
linkaxes([p1,p2],'xy');

%% 2.1
N = 1000;
V_Na = 55;
V_K  = -90;
n0 = 0.5;

V_vec = linspace(V_K+0.01, V_Na-0.01, N)';
n_nullcline = nullcline_n(V_vec);

v_nullcline = zeros(N,1);
for ii = 1:N
    v = V_vec(ii);
    v_nullcline(ii) = fzero(@(n)V_dot_2D(v, n), n0);
end

% plot
figure(4);
hold on
plot(V_vec, n_nullcline, 'r');
plot(V_vec, v_nullcline, 'b');
title("Nullclines for 2D HH approximation");
xlabel("V[mV]");
ylabel("n");
% plot([V_K, V_K], [0, max(v_nullcline)], 'g')
% plot([V_Na, V_Na], [0, max(v_nullcline)], 'g')
legend("$\dot{n}=0$", "$\dot{V}=0$", "interpreter", "latex") % "$V_K$", "$V_{Na}$"
ylim([0,1]);
hold off

%% 2.2
K = 100;
[N, V] = meshgrid(0:1/K:1, linspace(V_K+0.01, V_Na-0.01, K));
n_dots = n_dot(V, N);
V_dots = V_dot_2D(V, N);
n_plus = 2 * (n_dots > 0);
v_plus = V_dots > 0;

figure(6)
hold on
s = surf(V, N, n_plus + v_plus);
colormap lines;
view([0 90]);
xlabel("V[mV]");
ylabel("n");
xlim([V_K, V_Na]);
ylim([0,1]);
title("Signs of $\dot{V}$ and $\dot{n}$", "interpreter", "latex");
s.EdgeColor = 'none';
s1 = surf(V(1:2,1:2), N(1:2,1:2), ones(2));
s2 = surf(V(1:2,1:2), N(1:2,1:2), 2 * ones(2));
s3 = surf(V(1:2,1:2), N(1:2,1:2), zeros(2));
s4 = surf(V(1:2,1:2), N(1:2,1:2), 3 * ones(2));
legend( [s1, s2, s3, s4], ...
    "$\dot{n} < 0, \dot{V} < 0$", ...
    "$\dot{n} < 0, \dot{V} > 0$", ...
    "$\dot{n} > 0, \dot{V} < 0$", ...
    "$\dot{n} > 0, \dot{V} > 0$", ...
    "interpreter", "latex");

hold off;
%% 2.3
null_diff = (n_nullcline - v_nullcline) .^ 2;
zero_inds = find(null_diff < 1E-4);

ii = 1;
while ii <= length(zero_inds)
    if zero_inds(ii + 1) == zero_inds(ii) + 1
        if null_diff(zero_inds(ii)) < null_diff(zero_inds(ii + 1))
            zero_inds(ii+1) = [];
            ii = ii + 1;
        else
            zero_inds(ii) = [];
        end
    else
        ii = ii + 1;
    end 
end

V_min = V_vec(zero_inds);
n_min = v_nullcline(zero_inds);

delta = 1E-8;
n_dot_dot_V = (n_dot(V_min + delta, n_min) - n_dot(V_min - delta, n_min)) / (2 * delta);
V_dot_dot_V = (V_dot_2D(V_min + delta, n_min) - V_dot_2D(V_min - delta, n_min)) / (2 * delta);
n_dot_dot_n = (n_dot(V_min, n_min + delta) - n_dot(V_min, n_min - delta)) / (2 * delta);
V_dot_dot_n = (V_dot_2D(V_min, n_min + delta) - V_dot_2D(V_min, n_min - delta)) / (2 * delta);

J = zeros(2,2,3);
J(1,1,:) = V_dot_dot_V;
J(1,2,:) = V_dot_dot_n;
J(2,1,:) = n_dot_dot_V;
J(2,2,:) = n_dot_dot_n;

eig_vals = zeros(2,3);
for ii = 1:3
    eig_vals(:,ii) = eig(J(:,:,ii));
    point = [V_min(ii), n_min(ii)]
    lambda = eig_vals(:,ii)
end

%% 2.4
% use [t,y] = ode15s(@(t,y)hh_2D(t,y,0), tspan, y0);
n0 = 0.35 + x(7) / 100;
delta = 0.01;

v0 = -69;
v1 = -68;

N = 1E6;
n_vec1 = zeros(N, 1);
v_vec1 = zeros(N, 1);
n_vec1(1) = n0;
v_vec1(1) = v0;
n_vec2 = n_vec1;
v_vec2 = v_vec1;
v_vec2(1) = v1;
for ii = 2:N
    n_vec1(ii) = n_vec1(ii - 1) + n_dot(v_vec1(ii - 1), n_vec1(ii - 1)) * delta;
    v_vec1(ii) = v_vec1(ii - 1) + V_dot_2D(v_vec1(ii - 1), n_vec1(ii - 1)) * delta;
    n_vec2(ii) = n_vec2(ii - 1) + n_dot(v_vec2(ii - 1), n_vec2(ii - 1)) * delta;
    v_vec2(ii) = v_vec2(ii - 1) + V_dot_2D(v_vec2(ii - 1), n_vec2(ii - 1)) * delta;
end

figure(7);
hold on
plot(V_vec, n_nullcline);
plot(V_vec, v_nullcline);
plot(v_vec1, n_vec1, "LineWidth", 1.5);
plot(v_vec2, n_vec2, "LineWidth", 1.5);
title("Different starting points do or don't result in a PP");
xlabel("V[mV]");
ylabel("n");
legend("$\dot{n}=0$", "$\dot{V}=0$", ...
    "(V,n)=(" + num2str([v0,n0]) + ")", ...
    "(V,n)=(" + num2str([v1,n0]) + ")", ...
    "interpreter", "latex")
ylim([0,1]);
hold off
%% 2.5

%% 2.6

%% 2.7

%% 3.1

%% 3.2

%% 3.3

%% 3.4

%% 3.5

disp(">> done!");