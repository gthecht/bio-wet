close all; clear; clc;

%% Constants:
x = [2 0 6 8 3 7 0 0 7];
y = [2 0 3 8 2 7 2 9 0];
%% PART I

%% Definitions
gNa = 120; %[1/(K*Ohm*cm^2)]
gK = 36; %[1/(K*Ohm*cm^2)]
gCl = 0.3; %[1/(K*Ohm*cm^2)]

VNa = 55*10^-3; %[V]
VK = -90*10^-3; %[V]
VCl = -60*10^-3; %[V]

Cm = 1; %[uF/cm^2]

alpha_n = @(V_p) (0.1 - 0.01*((V_p + 74.4*10^-3)/(10^-3)))/(exp(1-0.1*((V_p + 74.4*10^-3)/(10^-3)))-1); %[1/ms]
alpha_m = @(V_p) (2.5 - 0.1*((V_p + 74.4*10^-3)/(10^-3)))/(exp(2.5-0.1*((V_p + 74.4*10^-3)/(10^-3)))-1); %[1/ms]
alpha_h = @(V_p) 0.07*exp(-((V_p + 74.4*10^-3)/(10^-3))/20); %[1/ms]
beta_n = @(V_p) 0.125*exp(-((V_p + 74.4*10^-3)/(10^-3))/80); %[1/ms]
beta_m = @(V_p) 4*exp(-((V_p + 74.4*10^-3)/(10^-3))/18); %[1/ms]
beta_h = @(V_p) 1/(exp(3-0.1*((V_p + 74.4*10^-3)/(10^-3)))+1); %[1/ms]

nInf = @(V_p) alpha_n(V_p)/(alpha_n(V_p) + beta_n(V_p));
mInf = @(V_p) alpha_m(V_p)/(alpha_m(V_p) + beta_m(V_p));
hInf = @(V_p) alpha_h(V_p)/(alpha_h(V_p) + beta_h(V_p));

I = 0; % Q1 assumption
VTag_rest = @(V_p) -(1/Cm)*(gNa*mInf(V_p)^3*hInf(V_p)*(V_p-VNa)+gK*nInf(V_p)^4*(V_p-VK)+gCl*(V_p-VCl))+I/Cm;


%% 1.1
V0 = fzero(VTag_rest, [VK VNa]);
n0 = nInf(V0);
m0 = mInf(V0);
h0 = hInf(V0);

disp(['1.1 the equilibrium point is: V = ' num2str(V0)]);
disp(['    n = ' num2str(n0) ', m = ' num2str(m0) ', h = ' num2str(h0)])

%% 1.2
syms n_p m_p h_p V_p
hh = [symfun(-(1/Cm)*(gNa*m_p^3*h_p*(V_p-VNa)+...
      gK*n_p^4*(V_p-VK)+gCl*(V_p-VCl))+I/Cm, [V_p, n_p, m_p, h_p]),...
      symfun(alpha_n*(1-n_p) - beta_n*n_p, [V_p, n_p, m_p, h_p]),...
      symfun(alpha_m*(1-m_p) - beta_m*m_p, [V_p, n_p, m_p, h_p]),...
      symfun(alpha_h*(1-h_p) - beta_h*h_p, [V_p, n_p, m_p, h_p])];

jacobian_sym = jacobian(hh, [V_p, n_p, m_p, h_p]);
jacobian_eig = eig(double(jacobian_sym(V0, n0, m0, h0)));

disp('1.2 The eigenvalues of the jacobian are: ');
disp(num2str(jacobian_eig));

%% 1.3
A = 1000;
y0 = [n0 m0 h0 V0];
tspan = 0:0.001:18;

T0 = 0.5697;
I0 = @(t) (10 + y(9))./(1+exp(A*(t - T0)))*10^-3;

%params = n m h V
hh0 = @(t, params) [alpha_n(params(4))*(1-params(1)) - beta_n(params(4))*params(1);...
                    alpha_m(params(4))*(1-params(2)) - beta_m(params(4))*params(2);...
                    alpha_h(params(4))*(1-params(3)) - beta_h(params(4))*params(3);...
                    -(1/Cm)*(gNa*params(2)^3*params(3)*(params(4)-VNa)+...
                            gK*params(1)^4*(params(4)-VK)+gCl*(params(4)-VCl))+I0(t)/Cm];

[t0,hhres0] = ode15s(hh0, tspan, y0);

T1 = T0+0.001;
I1 = @(t) (10 + y(9))./(1+exp(A*(t - T1)))*10^-3;

%params = n m h V
hh1 = @(t, params) [alpha_n(params(4))*(1-params(1)) - beta_n(params(4))*params(1);...
                    alpha_m(params(4))*(1-params(2)) - beta_m(params(4))*params(2);...
                    alpha_h(params(4))*(1-params(3)) - beta_h(params(4))*params(3);...
                    -(1/Cm)*(gNa*params(2)^3*params(3)*(params(4)-VNa)+...
                            gK*params(1)^4*(params(4)-VK)+gCl*(params(4)-VCl))+I1(t)/Cm];

[t1,hhres1] = ode15s(hh1, tspan, y0);

titles = {'n','m','h','V'};
yLab = {'n','m','h','V [Volt]'};
figure;
for index = 1:4
    subplot(2,2,index);
    hold on; grid on;
    plot(t0, hhres0(:,index), t1, hhres1(:,index));
    title(titles{index});
    legend('T0 = 0.5697ms' ,'T1 = 0.5707ms');
    xlabel('t[ms]');
    ylabel(yLab{index});
end

%% 2.1
I = 0;
n0 = 0.35 + y(8)/90;
h0 = n0;

VTag2D = @(V_p) -(1/Cm)*(gNa*mInf(V_p)^3*h0*(V_p-VNa)+gK*n0^4*(V_p-VK)+gCl*(V_p-VCl))+I/Cm;

Vdest = -(60 + y(7))*10^-3;
V02d = fzero(VTag2D, [Vdest-0.01 Vdest+0.01]);
m02d = mInf(V02d);

disp(['2.1 the equilibrium point is: V = ' num2str(V02d)]);
disp(['    n0 = h0 = ' num2str(n0) ', m0 = ' num2str(m02d)]);

syms m_p V_p
hh2 = [symfun(-(1/Cm)*(gNa*m_p^3*h0*(V_p-VNa)+gK*n0^4*(V_p-VK)+gCl*(V_p-VCl))+I/Cm, [V_p, m_p]),...
      symfun(alpha_m*(1-m_p) - beta_m*m_p, [V_p, m_p])];

jacobian_sym = jacobian(hh2, [V_p, m_p]);
jacobian_eig = eig(double(jacobian_sym(V02d, m02d)));

disp('2.2 The eigenvalues of the jacobian are:');
disp(num2str(jacobian_eig));

Vquiver = V02d - 2e-3 : 0.5e-3 : V02d + 2e-3;
mquiver = m02d - 0.035 : 0.005 : m02d + 0.035;

[VQ,MQ] = meshgrid(Vquiver,mquiver);
VRes = -(1/Cm)*(gNa.*MQ.^3*h0.*(VQ-VNa)+gK*n0^4.*(VQ-VK)+gCl.*(VQ-VCl));
mRes = (arrayfun(alpha_m,VQ).*(1-MQ) - arrayfun(beta_m,VQ).*MQ);

VRes_unit = (VRes.*1000)./sqrt((VRes.*1000).^2+(mRes.*100).^2);
mRes_unit = (mRes.*100)./sqrt((VRes.*1000).^2+(mRes.*100).^2);

figure;
scale = 0.23;
quiver(VQ*1000, MQ*100, VRes_unit, mRes_unit, scale);
title('State-Space Flow Lines Near Equilibrium Point of 2D HH Model');
xlabel('V [mV]');
ylabel('100*m');

%% 2.2
% No need for code

%% 3.1 + 3.2
t = [0.2 4];

for index = 1:2
    tspan = 0:0.001:t(index); %[ms]
    
    VstartPoint1 = V02d +2e-3;
    VstartPoint2 = V02d -2e-3;
    y0full1 = [n0 m02d h0 VstartPoint1];
    y0full2 = [n0 m02d h0 VstartPoint2];
    y02d1 = [m02d VstartPoint1];
    y02d2 = [m02d VstartPoint2];

    %params = n m h V
    hhFull = @(t, params) [alpha_n(params(4))*(1-params(1)) - beta_n(params(4))*params(1);...
                        alpha_m(params(4))*(1-params(2)) - beta_m(params(4))*params(2);...
                        alpha_h(params(4))*(1-params(3)) - beta_h(params(4))*params(3);...
                        -(1/Cm)*(gNa*params(2)^3*params(3)*(params(4)-VNa)+...
                                gK*params(1)^4*(params(4)-VK)+gCl*(params(4)-VCl))];

    [~,hhresFull1] = ode15s(hhFull, tspan, y0full1);
    [~,hhresFull2] = ode15s(hhFull, tspan, y0full2);

    % params = m V
    hh2d_sym = @(t,params) [alpha_m(params(2))*(1-params(1)) - beta_m(params(2))*params(1);...
        -(1/Cm)*(gNa*params(1)^3*h0*(params(2)-VNa)+ gK*n0^4*(params(2)-VK)+gCl*(params(2)-VCl))];

    [~,hhres2d1] = ode15s(hh2d_sym, tspan, y02d1);
    [~,hhres2d2] = ode15s(hh2d_sym, tspan, y02d2);

    figure;
    plot(V02d*1000,m02d,'*');
    hold on;
    grid on;
    plot(hhresFull1(:,4).*1000,hhresFull1(:,2), hhresFull2(:,4).*1000,hhresFull2(:,2),...
        hhres2d1(:,2).*1000, hhres2d1(:,1), hhres2d2(:,2).*1000, hhres2d2(:,1));
    xlabel('V [mV]');
    ylabel('m');
    title(['Simulation with t = ' num2str(t(index)) ' [ms]']);
    legend('Equilibrium Point','Full HH 1','Full HH 2', '2D HH 1', '2D HH 2');
    hold on;
end

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
figure;
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
figure;
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
K = 1000;
[N, V] = meshgrid(0:1/K:1, linspace(V_K+0.01, V_Na-0.01, K));
n_dots = n_dot(V, N);
V_dots = V_dot_2D(V, N);
n_plus = 2 * (n_dots > 0);
v_plus = V_dots > 0;

figure
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
s1 = surf(V(1:2,1:2), N(1:2,1:2), zeros(2));
s2 = surf(V(1:2,1:2), N(1:2,1:2), ones(2));
s3 = surf(V(1:2,1:2), N(1:2,1:2), 2 * ones(2));
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
[~, zero_inds] = findpeaks(-null_diff);

V_min = V_vec(zero_inds);
n_min = v_nullcline(zero_inds);

delta = 1E-8;
n_dot_d_V = (n_dot(V_min + delta, n_min) - n_dot(V_min - delta, n_min)) / (2 * delta);
V_dot_d_V = (V_dot_2D(V_min + delta, n_min) - V_dot_2D(V_min - delta, n_min)) / (2 * delta);
n_dot_d_n = (n_dot(V_min, n_min + delta) - n_dot(V_min, n_min - delta)) / (2 * delta);
V_dot_d_n = (V_dot_2D(V_min, n_min + delta) - V_dot_2D(V_min, n_min - delta)) / (2 * delta);

J = zeros(2,2,3);
J(1,1,:) = V_dot_d_V;
J(1,2,:) = V_dot_d_n;
J(2,1,:) = n_dot_d_V;
J(2,2,:) = n_dot_d_n;

eig_vals = zeros(2,3);
for ii = 1:3
    eig_vals(:,ii) = eig(J(:,:,ii));
    point = [V_min(ii), n_min(ii)]
    lambda = eig_vals(:,ii)
end

%% 2.4

n0 = 0.35 + x(7) / 100;
tspan = [0 500];
v0 = -69;
v1 = -68;

y0 = [n0, v0];
y1 = [n0, v1];

[t,y] = ode15s(@(t,y)hh2d(t,y,0), tspan, y0);
n_vec1 = y(:,1);
v_vec1 = y(:,2);

[t,y] = ode15s(@(t,y)hh2d(t,y,0), tspan, y1);
n_vec2 = y(:,1);
v_vec2 = y(:,2);

figure;
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

I = 3 + x(6) / 2;
N = 1000;
V_Na = 55;
V_K  = -90;
n0 = 0.5;

V_vec = linspace(V_K+0.01, V_Na-0.01, N)';
n_nullcline = nullcline_n(V_vec);

v_nullcline = zeros(N,1);
for ii = 1:N
    v = V_vec(ii);
    v_nullcline(ii) = fzero(@(n)V_dot_2D(v, n, I), n0);
end

% plot
figure;
hold on
plot(V_vec, n_nullcline, 'r');
plot(V_vec, v_nullcline, 'b');
title("Nullclines for 2D HH approximation");
xlabel("V[mV]");
ylabel("n");
% plot([V_K, V_K], [0, max(v_nullcline)], 'g')
% plot([V_Na, V_Na], [0, max(v_nullcline)], 'g')
legend("$\dot{n}=0$", "$\dot{V}=0$", "interpreter", "latex");
ylim([0,1]);

%% Now for the Jaacobians
null_diff = (n_nullcline - v_nullcline) .^ 2;
[~, zero_inds] = findpeaks(-null_diff);

V_min = V_vec(zero_inds);
n_min = v_nullcline(zero_inds);

text(V_min(1) - 3, n_min(1) + 0.03, "1");
text(V_min(2) - 3, n_min(2) + 0.03, "2");
text(V_min(3) - 3, n_min(3) + 0.03, "3");

delta = 1E-8;
n_dot_d_V = (n_dot(V_min + delta, n_min) - n_dot(V_min - delta, n_min)) / (2 * delta);
V_dot_d_V = (V_dot_2D(V_min + delta, n_min) - V_dot_2D(V_min - delta, n_min)) / (2 * delta);
n_dot_d_n = (n_dot(V_min, n_min + delta) - n_dot(V_min, n_min - delta)) / (2 * delta);
V_dot_d_n = (V_dot_2D(V_min, n_min + delta) - V_dot_2D(V_min, n_min - delta)) / (2 * delta);

J3 = zeros(2,2,3);
J3(1,1,:) = V_dot_d_V;
J3(1,2,:) = V_dot_d_n;
J3(2,1,:) = n_dot_d_V;
J3(2,2,:) = n_dot_d_n;

eig_vals3 = zeros(2,3);
for ii = 1:3
    eig_vals3(:,ii) = eig(J3(:,:,ii));
    point = [V_min(ii), n_min(ii)]
    lambda = eig_vals3(:,ii)
end

hold off

%% 3.2

%% 3.3
I_max = 8 - x(5) / 10;
I_min = 0;
k = 2;

while round(I_min, k) ~= round(I_max, k)
    I = (I_max + I_min) / 2;
    [points, lambdas] = get_equilibriums(I);
    ind = dsearchn(points', [-70, 0.35]);
    eig_vals = lambdas(:, ind);
    eig_val = real(eig_vals(1));
    if eig_val > 0
        I_max = I;
    elseif eig_val < 0
            I_min = I;
    end
end

I

%% 3.4
I = 3 + x(6) / 2;
N = 1000;
V_Na = 55;
V_K  = -90;
n0 = 0.5;

V_vec = linspace(V_K+0.01, V_Na-0.01, N)';
n_nullcline = nullcline_n(V_vec);

v_nullcline = zeros(N,1);
for ii = 1:N
    v = V_vec(ii);
    v_nullcline(ii) = fzero(@(n)V_dot_2D(v, n, I), n0);
end

K = 1000;
[N, V] = meshgrid(0:1/K:1.1, linspace(-100+0.01, 70-0.01, K));
n_dots = n_dot(V, N);
V_dots = V_dot_2D(V, N);
n_plus = 2 * (n_dots > 0);
v_plus = V_dots > 0;

figure
hold on
s = surf(V, N, n_plus + v_plus);
colormap lines;
view([0 90]);
xlabel("V[mV]");
ylabel("n");
xlim([-100, 70]);
ylim([0,1.1]);
title("Signs of $\dot{V}$ and $\dot{n}$", "interpreter", "latex");
s.EdgeColor = 'none';
s1 = surf(V(1:2,1:2), N(1:2,1:2), zeros(2));
s2 = surf(V(1:2,1:2), N(1:2,1:2), ones(2));
s3 = surf(V(1:2,1:2), N(1:2,1:2), 2 * ones(2));
s4 = surf(V(1:2,1:2), N(1:2,1:2), 3 * ones(2));
p1 = plot3([-90,-90,60,60,-90],[0,1,1,0,0], [8,8,8,8,8], 'g', "LineWidth", 1.5);
legend( [s1, s2, s3, s4, p1], ...
    "$\dot{n} < 0, \dot{V} < 0$", ...
    "$\dot{n} < 0, \dot{V} > 0$", ...
    "$\dot{n} > 0, \dot{V} < 0$", ...
    "$\dot{n} > 0, \dot{V} > 0$", ...
    "Outer Ring", ...
    "interpreter", "latex");
hold off;

%% 3.5
I = 5 + x(4);
N = 1000;
V_Na = 55;
V_K  = -90;
n0 = 0.5;

V_vec = linspace(V_K+0.01, V_Na-0.01, N)';
n_nullcline = nullcline_n(V_vec);

v_nullcline = zeros(N,1);
for ii = 1:N
    v = V_vec(ii);
    v_nullcline(ii) = fzero(@(n)V_dot_2D(v, n, I), n0);
end

y0 = [0.25, -70];
y1 = [0.45, -70];

tspan = [0 500];

[t,y] = ode15s(@(t,y)hh2d(t,y,I), tspan, y0);
n_vec1 = y(:,1);
v_vec1 = y(:,2);

[t,y] = ode15s(@(t,y)hh2d(t,y,I), tspan, y1);
n_vec2 = y(:,1);
v_vec2 = y(:,2);


figure;
hold on
plot(V_vec, n_nullcline);
plot(V_vec, v_nullcline);
plot(v_vec1, n_vec1, "LineWidth", 1.5);
plot(v_vec2, n_vec2, "LineWidth", 1.5);
title("Different starting points do or don't result in a PP");
xlabel("V[mV]");
ylabel("n");
legend("$\dot{n}=0$", "$\dot{V}=0$", ...
    "(V,n)=(" + num2str([y0(1),y0(2)]) + ")", ...
    "(V,n)=(" + num2str([y1(1), y1(2)]) + ")", ...
    "interpreter", "latex")
ylim([0,1]);
hold off

disp(">> done!");