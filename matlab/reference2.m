clear; close all;

%% Wet Assignment- Hodgkin - Huxley

% Submitters:
% Yael David
x = [2 0 4 8 6 1 4 1 3];
% Shunit Haviv
y = [3 0 5 1 4 7 4 5 6];

%% Definitions
gNa = 120; %[1/(K*Ohm*cm^2)]
gK = 36; %[1/(K*Ohm*cm^2)]
gCl = 0.3; %[1/(K*Ohm*cm^2)]

VNa = 55*10^-3; %[V]
VK = -90*10^-3; %[V]
VCl = -60*10^-3; %[V]

Cm = 1; %[uF/cm^2]

alpha_n = @(V) (0.1 - 0.01*((V + 74.4*10^-3)/(10^-3)))/(exp(1-0.1*((V + 74.4*10^-3)/(10^-3)))-1); %[1/ms]
alpha_m = @(V) (2.5 - 0.1*((V + 74.4*10^-3)/(10^-3)))/(exp(2.5-0.1*((V + 74.4*10^-3)/(10^-3)))-1); %[1/ms]
alpha_h = @(V) 0.07*exp(-((V + 74.4*10^-3)/(10^-3))/20); %[1/ms]
beta_n = @(V) 0.125*exp(-((V + 74.4*10^-3)/(10^-3))/80); %[1/ms]
beta_m = @(V) 4*exp(-((V + 74.4*10^-3)/(10^-3))/18); %[1/ms]
beta_h = @(V) 1/(exp(3-0.1*((V + 74.4*10^-3)/(10^-3)))+1); %[1/ms]

tau_n = @(V) 1/(alpha_n(V) + beta_n(V));
tau_m = @(V) 1/(alpha_m(V) + beta_m(V));
tau_h = @(V) 1/(alpha_h(V) + beta_h(V));
nInf = @(V) alpha_n(V)/(alpha_n(V) + beta_n(V));
mInf = @(V) alpha_m(V)/(alpha_m(V) + beta_m(V));
hInf = @(V) alpha_h(V)/(alpha_h(V) + beta_h(V));

I = 0; %assumption Q1
VTag_rest = @(V) -(1/Cm)*(gNa*mInf(V)^3*hInf(V)*(V-VNa)+gK*nInf(V)^4*(V-VK)+gCl*(V-VCl))+I/Cm;

%% Part A
disp('Part A');
%% Q1
%% Q1.1

V0 = fzero(VTag_rest, [VK VNa]);
n0 = nInf(V0);
m0 = mInf(V0);
h0 = hInf(V0);

disp(['1.1 the equilibrium point is: V = ' num2str(V0)]);
disp(['    n = ' num2str(n0) ', m = ' num2str(m0) ', h = ' num2str(h0)])

%% Q1.2

syms n m h V
hh = [symfun(-(1/Cm)*(gNa*m^3*h*(V-VNa)+...
      gK*n^4*(V-VK)+gCl*(V-VCl))+I/Cm, [V, n, m, h]),...
      symfun(alpha_n*(1-n) - beta_n*n, [V, n, m, h]),...
      symfun(alpha_m*(1-m) - beta_m*m, [V, n, m, h]),...
      symfun(alpha_h*(1-h) - beta_h*h, [V, n, m, h])];

jacobian_sym = jacobian(hh, [V, n, m, h]);
jacobian_eig = eig(double(jacobian_sym(V0, n0, m0, h0)));

disp('1.2 The eigenvalues of the jacobian are:');
disp(num2str(jacobian_eig));

%% Q1.3

A = 1000;
y0 = [n0 m0 h0 V0];
tspan = 0:0.001:18;

T0 = 0.3995;
I0 = @(t) (10 + y(9))./(1+exp(A*(t - T0)))*10^-3;

%params = n m h V
hh0 = @(t, params) [alpha_n(params(4))*(1-params(1)) - beta_n(params(4))*params(1);...
                    alpha_m(params(4))*(1-params(2)) - beta_m(params(4))*params(2);...
                    alpha_h(params(4))*(1-params(3)) - beta_h(params(4))*params(3);...
                    -(1/Cm)*(gNa*params(2)^3*params(3)*(params(4)-VNa)+...
                            gK*params(1)^4*(params(4)-VK)+gCl*(params(4)-VCl))+I0(t)/Cm];

[t0,hhres0] = ode15s(hh0, tspan, y0);

T1 = 0.3996;
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
    plot(t0, hhres0(:,index), t1, hhres1(:,index));
    title(titles{index});
    legend('T = 0.3995ms','T = 0.3996ms');
    xlabel('t[ms]');
    ylabel(yLab{index});
end

%% Q2
%% Q2.1
n0 = 0.35 + y(8)/90;
h0 = n0;

VTag2D = @(V) -(1/Cm)*(gNa*mInf(V)^3*h0*(V-VNa)+gK*n0^4*(V-VK)+gCl*(V-VCl))+I/Cm;

Vdest = -(60 + y(7))*10^-3;
V02d = fzero(VTag2D, [Vdest-0.01 Vdest+0.01]);
m02d = mInf(V02d);

disp(['2.1 the equilibrium point is: V = ' num2str(V02d)]);
disp(['    n0 = h0 = ' num2str(n0) ', m0 = ' num2str(m02d)]);

syms m V
hh2 = [symfun(-(1/Cm)*(gNa*m^3*h0*(V-VNa)+gK*n0^4*(V-VK)+gCl*(V-VCl))+I/Cm, [V, m]),...
      symfun(alpha_m*(1-m) - beta_m*m, [V, m])];

jacobian_sym = jacobian(hh2, [V, m]);
jacobian_eig = eig(double(jacobian_sym(V02d, m02d)));

disp('2.2 The eigenvalues of the jacobian are:');
disp(num2str(jacobian_eig));

Vquiver = V02d - 2e-3 : 0.5e-3 : V02d + 2e-3;
mquiver = m02d - 0.035 : 0.005 : m02d + 0.035;

[MQ, VQ] = meshgrid(mquiver,Vquiver);
VRes = -(1/Cm)*(gNa.*MQ.^3*h0.*(VQ-VNa)+gK*n0^4.*(VQ-VK)+gCl.*(VQ-VCl));
mRes = (arrayfun(alpha_m,VQ).*(1-MQ) - arrayfun(beta_m,VQ).*MQ);

figure;
quiver(MQ, VQ, mRes, VRes);
xlabel('m');
ylabel('V [Volt]');

%% Q2.2, 3.1
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
    hh2d = @(t,params) [alpha_m(params(2))*(1-params(1)) - beta_m(params(2))*params(1);...
        -(1/Cm)*(gNa*params(1)^3*h0*(params(2)-VNa)+ gK*n0^4*(params(2)-VK)+gCl*(params(2)-VCl))];

    [~,hhres2d1] = ode15s(hh2d, tspan, y02d1);
    [~,hhres2d2] = ode15s(hh2d, tspan, y02d2);

    figure;
    plot(V02d,m02d,'*');
    hold on;
    plot(hhresFull1(:,4),hhresFull1(:,2), hhresFull2(:,4),hhresFull2(:,2),...
        hhres2d1(:,2), hhres2d1(:,1), hhres2d2(:,2), hhres2d2(:,1));
    xlabel('V [Volt]');
    ylabel('m');
    title(['simulation with t = ' num2str(t(index)) 'ms']);
    legend('Equilibrium point','Full HH 1','Full HH 2', '2D HH 1', '2D HH 2');
    hold on;

end

%% Part B
disp('Part B');
close all;
%% Q1

m0 = 0.05;
n0 = 0.32;
h0 = 0.6;
V0 = -74*10^-3;
y0 = [n0 m0 h0 V0];
tspan = 0:0.05:200;

% a = [q1.1 q1.2 q1.4 q1.4]
a = [(1+x(9)/9) (5+x(8)) 4 2.271 4.2];

for index = 1:numel(a)

    I0 = @(t) a(index)*10^-3;

    %params = n m h V
    hh3 = @(t, params) [alpha_n(params(4))*(1-params(1)) - beta_n(params(4))*params(1);...
                        alpha_m(params(4))*(1-params(2)) - beta_m(params(4))*params(2);...
                        alpha_h(params(4))*(1-params(3)) - beta_h(params(4))*params(3);...
                        -(1/Cm)*(gNa*params(2)^3*params(3)*(params(4)-VNa)+...
                                gK*params(1)^4*(params(4)-VK)+gCl*(params(4)-VCl))+I0(t)/Cm];

    [t0,hhres0] = ode15s(hh3, tspan, y0);

    titles = {'n','m','h','V'};
    yLab = {'n','m','h','V [Volt]'};
    figure;
    for ind = 1:4
        subplot(2,2,ind);
        plot(t0, hhres0(:,ind));
        title([titles{ind} ' with a = ' num2str(a(index)) '\muA']);
        xlabel('t[ms]');
        ylabel(yLab{ind});
    end
end

%% Q2
%% Q2.1

vRange = VK+0.001:0.001:VNa-0.001;

I = 0;
VTag3 = symfun(-(1/Cm)*(gNa*mInf(V)^3*(1-n)*(V-VNa)+gK*n^4*(V-VK)+gCl*(V-VCl))+I/Cm, [V, n]);
nTag3 = symfun(alpha_n(V)*(1-n) - beta_n(V)*n, [V, n]);

nNull = solve(nTag3(V,n) == 0,n);
nNullcline = double(subs(nNull,V,vRange));

vNull = solve(VTag3(V,n),n);
vNull = vNull(2);

for ind = 1:numel(vRange)
    vNullcline(ind) = double(subs(vNull,V,vpa(vRange(ind))));
end

figure;
plot(vRange,nNullcline,vRange,vNullcline);
title('V,n Nullclines');
xlabel('V [Volt]');
ylabel('n');
legend('n Nullcline','V Nullcline');

%% Q2.3

clear VTagVec;
for ind = 1:numel(vRange)
    VTagVec(ind) = vpa(subs(VTag3,{V, n}, {vRange(ind), nNullcline(ind)}));
end

VTagVec(VTagVec<0) = 0;
VTagVec(VTagVec>0) = 1;

indexList = [];
for ind = 1:numel(VTagVec) -1
    if VTagVec(ind) ~= VTagVec(ind + 1)
        indexList = [indexList ind];
    end
end

vEq = vRange(indexList);
nEq = nNullcline(indexList);

syms m V
hh4 = [symfun(-(1/Cm)*(gNa*mInf(V)^3*(1-n)*(V-VNa)+gK*n^4*(V-VK)+gCl*(V-VCl))+I/Cm, [V, n]),...
      symfun(alpha_n(V)*(1-n) - beta_n(V)*n, [V, n])];  
  
jacobian_sym3 = jacobian(hh4, [V, n]);

for ind = 1:3
    jacobian_eig3(ind,:) = eig(double(jacobian_sym3(vEq(ind), nEq(ind))));
end

disp('2.3 The equilibeium points and the eigenvalues of the jacobian are:');
disp('  V           n                       Eigenvalues');
spaces = ['     '; '     '; '     '];
disp([num2str(vEq') spaces num2str(nEq') spaces num2str(jacobian_eig3)]);


%% Q2.4
close all;

figure;
plot(vRange,nNullcline,vRange,vNullcline);
xlabel('V [Volt]');
ylabel('n');
hold on;


nInit = 0.35 + x(7)/100;
VInit = [-0.0665 -0.067];

tspan = 0:0.01:200;

for index = 1:numel(VInit)

    y0 = [VInit(index) nInit];  
    
    %params = V n
    hh5 = @(t, params) [-(1/Cm)*(gNa*mInf(params(1))^3*(1-params(2))*(params(1)-VNa)+...
                        gK*params(2)^4*(params(1)-VK)+gCl*(params(1)-VCl));...
                        alpha_n(params(1))*(1-params(2)) - beta_n(params(1))*params(2)];

    [~,hhres3] = ode15s(hh5, tspan, y0);

    plot(hhres3(:,1),hhres3(:,2));
end
legend('n Nullcline','V Nullcline','Action potential', 'No action potential');

%% Q3
%% Q3.1
close all;

vRange = VK+0.001:0.001:VNa-0.001;
nRange = 0.001:0.005:0.99;

I = (3 + x(6)/2)*10^-3; % I = 0.0035;
% I = 3.33*10^-3;  %Bifurcation value
VTag3 = symfun(-(1/Cm)*(gNa*mInf(V)^3*(1-n)*(V-VNa)+gK*n^4*(V-VK)+gCl*(V-VCl))+I/Cm, [V, n]);
nTag3 = symfun(alpha_n(V)*(1-n) - beta_n(V)*n, [V, n]);

nNull = solve(isolate(nTag3(V,n) == 0,n),n);
nNullcline = double(subs(nNull,V,vRange));

vNull = solve(VTag3(V,n),n);
vNull = vNull(2);

for ind = 1:numel(vRange)
    vNullcline(ind) = double(subs(vNull,V,vpa(vRange(ind))));
end

figure;
plot(vRange,nNullcline,vRange,vNullcline);
title('V,n Nullclines');
xlabel('V [Volt]');
ylabel('n');
legend('n Nullcline','V Nullcline');

%%Q3.1

clear VTagVec;
for ind = 1:numel(vRange)
    VTagVec(ind) = vpa(subs(VTag3,{V, n}, {vRange(ind), nNullcline(ind)}));
end

VTagVec(VTagVec<0) = 0;
VTagVec(VTagVec>0) = 1;

indexList = [];
for ind = 1:numel(VTagVec) -1
    if VTagVec(ind) ~= VTagVec(ind + 1)
        indexList = [indexList ind];
    end
end

vEq = vRange(indexList);
nEq = nNullcline(indexList);

syms m V
hh4 = [symfun(-(1/Cm)*(gNa*mInf(V)^3*(1-n)*(V-VNa)+gK*n^4*(V-VK)+gCl*(V-VCl))+I/Cm, [V, n]),...
      symfun(alpha_n(V)*(1-n) - beta_n(V)*n, [V, n])];  
  
jacobian_sym3 = jacobian(hh4, [V, n]);

for ind = 1:3
    jacobian_eig3(ind,:) = eig(double(jacobian_sym3(vEq(ind), nEq(ind))));
end

disp('2.3 The equilibeium points and the eigenvalues of the jacobian are:');
disp('  V           n                       Eigenvalues');
spaces = ['     '; '     '; '     '];
disp([num2str(vEq') spaces num2str(nEq') spaces num2str(jacobian_eig3)]);


%% Q3.5
close all;

figure;
plot(vRange,nNullcline,vRange,vNullcline);
xlabel('V [Volt]');
ylabel('n');
hold on;

nInit = [0.25 0.45];
VInit = -70*10^-3;

tspan = 0:0.01:200;

for index = 1:numel(nInit)

    y0 = [VInit nInit(index)];  
    
    %params = V n
    hh5 = @(t, params) [-(1/Cm)*(gNa*mInf(params(1))^3*(1-params(2))*(params(1)-VNa)+...
                        gK*params(2)^4*(params(1)-VK)+gCl*(params(1)-VCl)) + I/Cm;...
                        alpha_n(params(1))*(1-params(2)) - beta_n(params(1))*params(2)];

    [~,hhres3] = ode15s(hh5, tspan, y0);

    plot(hhres3(:,1),hhres3(:,2));
end

legend('n Nullcline','V Nullcline','n = 25', 'n = 45');
