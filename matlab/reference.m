%% Tal Bachar 208115485 Yotam David 206718777
%% part a
clear all
close all
%% 1.1
gna = 120;
gk = 36;
gcl = 0.3;
vna = 55;
vk = -90;
vcl = -60;
cm = 1;
u = @(v) v+74.44;
a_n = @(u) (0.1-0.01*u+eps)/(exp(1-0.1*u)-1+eps);
a_m = @(u) (2.5-0.1*u+eps)/(exp(2.5-0.1*u)-1+eps);
a_h = @(u) 0.07*exp(-u/20)+eps;
b_n = @(u) 0.125*exp(-u/80)+eps;
b_m = @(u) 4*exp(-u/18)+eps;
b_h = @(u) (1+eps)/(exp(3-0.1*u)+1+eps);

n_inf = @(u) a_n(u)/(a_n(u)+b_n(u));
m_inf = @(u) a_m(u)/(a_m(u)+b_m(u));
h_inf = @(u) a_h(u)/(a_h(u)+b_h(u));
vdot_inf = @(u) -(1/cm)*(gna*(m_inf(u)^3)*h_inf(u)*(u-74.44-vna)+gk*(n_inf(u)^4)*(u-74.44-vk)+gcl*(u-74.44-vcl));

u_zero = fzero(vdot_inf,0);
equalibrium = [m_inf(u_zero),n_inf(u_zero),h_inf(u_zero),u_zero];
%% 1.2
%phase state equations
syms m n h u
mdot(m,n,h,u) =  a_m(u)*(1-m)-b_m(u)*m; % dm/dt
ndot(m,n,h,u) =  a_n(u)*(1-n)-b_n(u)*n; % dn/dt
hdot(m,n,h,u) =  a_h(u)*(1-h)-b_h(u)*h; % dh/dt
vdot(m,n,h,u) =  -(1/cm)*(gna*(m^3)*h*(u-74.44-vna)+gk*(n^4)*(u-74.44-vk)+gcl*(u-74.44-vcl)); % dv/dt
%Jacobian
symvec = [m,n,h,u];
functionvec = [mdot(m,n,h,u),ndot(m,n,h,u),hdot(m,n,h,u),vdot(m,n,h,u)];
J = zeros(4,4);
for ii=1:4
    for jj=1:4
        df = symfun(diff(functionvec(ii),symvec(jj)),symvec);
        J(ii,jj) = double(df(equalibrium(1),equalibrium(2),equalibrium(3),equalibrium(4)));
    end
end
lambda = eig(J);
%% 1.3
T=0.382;
t = 0:0.1:18;
I= @(t) 15./(1+exp(1000*(t-T)))-15./(1+exp(1000*(t)));
y0 = [equalibrium(1),equalibrium(2),equalibrium(3),equalibrium(4)-74.44];
[t,y] = ode15s(@(t,y)hhx(t,y,I),t,y0);

figure(1);
subplot(2,2,1);
plot(t,y(:,1));
ylabel("m");
xlabel("t[ms]");
hold on;

subplot(2,2,2);
plot(t,y(:,2));
ylabel("n");
xlabel("t[ms]");
hold on;

subplot(2,2,3);
plot(t,y(:,3));
ylabel("h");
xlabel("t[ms]");
hold on;

subplot(2,2,4);
plot(t,y(:,4));
ylabel("V[mv]");
xlabel("t[ms]");
hold on;

T=0.381;
I= @(t) 15./(1+exp(1000*(t-T)))-15./(1+exp(1000*(t)));
y0 = [equalibrium(1),equalibrium(2),equalibrium(3),equalibrium(4)-74.44];
[t,y] = ode15s(@(t,y)hhx(t,y,I),t,y0);

subplot(2,2,1);
plot(t,y(:,1));
hold off;

subplot(2,2,2);
plot(t,y(:,2));
hold off;

subplot(2,2,3);
plot(t,y(:,3));
hold off;

subplot(2,2,4);
plot(t,y(:,4));
hold off;

%% 2.1
h0=0.35+8/90;
n0=h0;
r0=-64+74.44;
vdot_inf_2 = @(u) -(1/cm)*(gna*(m_inf(u)^3)*h0*(u-74.44-vna)+gk*(n0^4)*(u-74.44-vk)+gcl*(u-74.44-vcl));
u_zero_2 = fzero(vdot_inf_2,r0);
equalibrium_2 = [m_inf(u_zero_2),n0,h0,u_zero_2];

syms m u
mdot_2(m,u) =  a_m(u)*(1-m)-b_m(u)*m; % dm/dt
vdot_2(m,u) =  -(1/cm)*(gna*(m_inf(u)^3)*h0*(u-74.44-vna)+gk*(n0^4)*(u-74.44-vk)+gcl*(u-74.44-vcl)); % dv/dt
%Jacobian
symvec_2 = [m,u];
functionvec_2 = [mdot_2(m,u),vdot_2(m,u)];
J2 = zeros(2,2);
for ii=1:2
    for jj=1:2
        df = symfun(diff(functionvec_2(ii),symvec_2(jj)),symvec_2);
        J2(ii,jj) = double(df(equalibrium_2(1),equalibrium_2(4)));
    end
end
lambda2 = eig(J2);
%quiver
around_u=(equalibrium_2(4)-4):0.5:(equalibrium_2(4)+4);
around_m0=(equalibrium_2(1)-0.07):0.005:(equalibrium_2(1)+0.07);
[u1,m1]= meshgrid(around_u,around_m0);
am = (2.5-0.1*u1+eps)./(exp(2.5-0.1*u1)-1+eps);
bm = 4*exp(-u1/18)+eps;
dm1 =  am.*(1-m1)-bm.*m1; % dm/dt
dv1 =  -(1/cm)*(gna*(m1.^3).*h0.*(u1-74.44-vna)+gk.*(n0.^4).*(u1-74.44-vk)+gcl.*(u1-74.44-vcl)); % dv/dt
figure(2);
quiver(u1-74.44,100*m1,dv1,100*dm1);
title('equlibrium center in (v0,m0) field lines');
xlabel('v');
ylabel('100m');
hold on;
%% 3.1
dt=0.001;
maxT=0.2;
t=0:dt:maxT;
init=[equalibrium_2(4)+2,equalibrium_2(4)-2];

for ii=1:2
x=zeros(1,maxT/dt+1);
x(1)=init(ii);
y=zeros(1,maxT/dt+1);
y(1)=equalibrium_2(1);

f = @(x,y)  -(1/cm)*(gna*(y^3)*h0*(x-74.44-vna)+gk*(n0^4)*(x-74.44-vk)+gcl*(x-74.44-vcl));
g = @(x,y)  a_m(x)*(1-y)-b_m(x)*y;

for k=1:maxT/dt
    x(k+1)=x(k)+f(x(k),y(k))*dt;
    y(k+1)=y(k)+g(x(k),y(k))*dt;
end
figure(3);
p(1)=plot(x-74.44,y,'color','blue');
hold on;
end
%hold off;

%full HH
init=[equalibrium_2(4)+2,equalibrium_2(4)-2];

for ii=1:2
x=zeros(1,maxT/dt+1);
x(1)=init(ii);
y=zeros(1,maxT/dt+1);
y(1)=equalibrium_2(1);
z=zeros(1,maxT/dt+1);
z(1)=equalibrium_2(2);
w=zeros(1,maxT/dt+1);
w(1)=equalibrium_2(3);

f = @(x,y,z,w)  -(1/cm)*(gna*(y^3)*w*(x-74.44-vna)+gk*(z^4)*(x-74.44-vk)+gcl*(x-74.44-vcl));
s = @(x,z)  a_n(x)*(1-z)-b_n(x)*z;
l = @(x,w)  a_h(x)*(1-w)-b_h(x)*w;
g = @(x,y)  a_m(x)*(1-y)-b_m(x)*y;

for k=1:maxT/dt
    x(k+1)=x(k)+f(x(k),y(k),z(k),w(k))*dt;
    y(k+1)=y(k)+g(x(k),y(k))*dt;
    z(k+1)=z(k)+s(x(k),z(k))*dt;
    w(k+1)=w(k)+l(x(k),w(k))*dt;
end
%figure(4);
p(2) = plot(x-74.44,y,'color','red');
hold on;
end
hold off;
xlabel("V[mv]");
ylabel("m");
legend([p(1) p(2)],"2D model","Full model");
%%
dt=0.001;
maxT=4;
t=0:dt:maxT;
init=[equalibrium_2(4)+2,equalibrium_2(4)-2];

for ii=1:2
x=zeros(1,maxT/dt+1);
x(1)=init(ii);
y=zeros(1,maxT/dt+1);
y(1)=equalibrium_2(1);

f = @(x,y)  -(1/cm)*(gna*(y^3)*h0*(x-74.44-vna)+gk*(n0^4)*(x-74.44-vk)+gcl*(x-74.44-vcl));
g = @(x,y)  a_m(x)*(1-y)-b_m(x)*y;

for k=1:maxT/dt
    x(k+1)=x(k)+f(x(k),y(k))*dt;
    y(k+1)=y(k)+g(x(k),y(k))*dt;
end
figure(4);
p(1)=plot(x-74.44,y,'color','blue');
hold on;
end
%hold off;

%full HH
init=[equalibrium_2(4)+2,equalibrium_2(4)-2];

for ii=1:2
x=zeros(1,maxT/dt+1);
x(1)=init(ii);
y=zeros(1,maxT/dt+1);
y(1)=equalibrium_2(1);
z=zeros(1,maxT/dt+1);
z(1)=equalibrium_2(2);
w=zeros(1,maxT/dt+1);
w(1)=equalibrium_2(3);

f = @(x,y,z,w)  -(1/cm)*(gna*(y^3)*w*(x-74.44-vna)+gk*(z^4)*(x-74.44-vk)+gcl*(x-74.44-vcl));
s = @(x,z)  a_n(x)*(1-z)-b_n(x)*z;
l = @(x,w)  a_h(x)*(1-w)-b_h(x)*w;
g = @(x,y)  a_m(x)*(1-y)-b_m(x)*y;

for k=1:maxT/dt
    x(k+1)=x(k)+f(x(k),y(k),z(k),w(k))*dt;
    y(k+1)=y(k)+g(x(k),y(k))*dt;
    z(k+1)=z(k)+s(x(k),z(k))*dt;
    w(k+1)=w(k)+l(x(k),w(k))*dt;
end
%figure(4);
p(2) = plot(x-74.44,y,'color','red');
hold on;
end
hold off;
xlabel("V[mv]");
ylabel("m");
legend([p(1) p(2)],"2D model","Full model");
%% part b
%% 1.1
t = 0:0.1:200;
I= @(t) (1+7/9)*(1-(1./(1+exp(1000*(t)))));
y0 = [equalibrium(1),equalibrium(2),equalibrium(3),equalibrium(4)-74.44];
[t,y] = ode15s(@(t,y)hhx(t,y,I),t,y0);

figure(5);
subplot(2,2,1);
plot(t,y(:,1));
ylabel("m");
xlabel("t[ms]");
hold on;

subplot(2,2,2);
plot(t,y(:,2));
ylabel("n");
xlabel("t[ms]");
hold on;

subplot(2,2,3);
plot(t,y(:,3));
ylabel("h");
xlabel("t[ms]");
hold on;

subplot(2,2,4);
plot(t,y(:,4));
ylabel("V[mv]");
xlabel("t[ms]");
hold off;
%% 1.2
t = 0:0.1:200;
I= @(t) (5+7)*(1-(1./(1+exp(1000*(t)))));
y0 = [equalibrium(1),equalibrium(2),equalibrium(3),equalibrium(4)-74.44];
[t,y] = ode15s(@(t,y)hhx(t,y,I),t,y0);

figure(6);
subplot(2,2,1);
plot(t,y(:,1));
ylabel("m");
xlabel("t[ms]");
hold on;

subplot(2,2,2);
plot(t,y(:,2));
ylabel("n");
xlabel("t[ms]");
hold on;

subplot(2,2,3);
plot(t,y(:,3));
ylabel("h");
xlabel("t[ms]");
hold on;

subplot(2,2,4);
plot(t,y(:,4));
ylabel("V[mv]");
xlabel("t[ms]");
hold off;
%% 1.4
t = 0:0.1:200;
I= @(t) 4.1*(1-(1./(1+exp(1000*(t)))));
y0 = [equalibrium(1),equalibrium(2),equalibrium(3),equalibrium(4)-74.44];
[t,y] = ode15s(@(t,y)hhx(t,y,I),t,y0);

figure(7);
subplot(2,2,1);
plot(t,y(:,1));
ylabel("m");
xlabel("t[ms]");
hold on;

subplot(2,2,2);
plot(t,y(:,2));
ylabel("n");
xlabel("t[ms]");
hold on;

subplot(2,2,3);
plot(t,y(:,3));
ylabel("h");
xlabel("t[ms]");
hold on;

subplot(2,2,4);
plot(t,y(:,4));
ylabel("V[mv]");
xlabel("t[ms]");
hold off;
%% 2.1
n_inf_values = zeros(1,vna-vk+1);
for ii=vk:vna
n_inf_values(ii-vk+1) = n_inf(ii+74.44);
end
figure(8);
title("phase plane nullclines");
xlabel("V[mv]");
ylabel("n");
plot(vk:1:vna,n_inf_values,'--','color','blue');
hold on;
v_inf_values = zeros(1,vna-vk+1);
for ii=vk+0.25:1:vna
v_dot_inf = @(n) -(1/cm)*(gna*(m_inf(ii+74.44)^3)*(1-n)*(ii-vna)+gk*(n^4)*(ii-vk)+gcl*(ii-vcl));
v_inf_values(ii-(vk+0.25)+1) = fzero(v_dot_inf,[0 1]);
end
plot(vk:1:vna,v_inf_values(1:146),'-.','color','red');

hold on;
%% 2.3
n_zero_points = [17,32,63];
v_zero_points = [-74,-59,-29];
syms n u
n_dot(n,u) =  a_n(u)*(1-n)-b_n(u)*n; % dm/dt
v_dot(n,u) =  -(1/cm)*(gna*(m_inf(u)^3)*(1-n)*(u-74.44-vna)+gk*(n^4)*(u-74.44-vk)+gcl*(u-74.44-vcl)); % dv/dt
%Jacobian
symvec_3 = [n,u];
functionvec_3 = [n_dot(n,u),v_dot(n,u)];
lambda3 = zeros(3,2);
for kk=1:3
J3 = zeros(2,2);
for ii=1:2
    for jj=1:2
        df = symfun(diff(functionvec_3(ii),symvec_3(jj)),symvec_3);
        J3(ii,jj) = double(df(n_inf_values(n_zero_points(kk)),v_zero_points(kk)+74.44));
    end
end
lambda3(kk,:) = eig(J3);
end
%% 2.4
t = 0:0.1:100;
n_0 = 0.42;
v_0 = -67;
I= @(t) 0;
y0 = [n_0,v_0];
[t,y] = ode15s(@(t,y)hhx2(t,y,I),t,y0);

plot(y(:,2),y(:,1),'color','magenta');
hold on;

v_0 = -65;
y0 = [n_0,v_0];
[t,y] = ode15s(@(t,y)hhx2(t,y,I),t,y0);

plot(y(:,2),y(:,1),'color','green');
legend('n-nullcline','v-nullcline','v0=-67','v0=-65');
hold off;
%% 3.1
I_0 = 12;
n_inf_values = zeros(1,vna-vk+1);
for ii=vk:vna
n_inf_values(ii-vk+1) = n_inf(ii+74.44);
end
figure(9);
title("solution courses in the phase plane");
xlabel("V[mv]");
ylabel("n");
plot(vk:1:vna,n_inf_values,'--','color','blue');
hold on;
v_inf_values = zeros(1,vna-vk+1);
for ii=vk+0.6:1:vna
v_dot_inf = @(n) -(1/cm)*(gna*(m_inf(ii+74.44)^3)*(1-n)*(ii-vna)+gk*(n^4)*(ii-vk)+gcl*(ii-vcl))+I_0/cm;
v_inf_values(ii-(vk+0.6)+1) = fzero(v_dot_inf,1);
end
plot(vk:1:vna,v_inf_values(1:146),'-.','color','red');
hold on;
n_zero_points = [20,29,62];
v_zero_points = [-71,-61,-28];
syms n u
n_dot(n,u) =  a_n(u)*(1-n)-b_n(u)*n; % dm/dt
v_dot(n,u) =  -(1/cm)*(gna*(m_inf(u)^3)*(1-n)*(u-74.44-vna)+gk*(n^4)*(u-74.44-vk)+gcl*(u-74.44-vcl))+I_0/cm; % dv/dt
%Jacobian
symvec_4 = [n,u];
functionvec_4 = [n_dot(n,u),v_dot(n,u)];
lambda4 = zeros(3,2);
for kk=1:3
J4 = zeros(2,2);
for ii=1:2
    for jj=1:2
        df = symfun(diff(functionvec_4(ii),symvec_4(jj)),symvec_4);
        J4(ii,jj) = double(df(n_inf_values(n_zero_points(kk)),v_zero_points(kk)+74.44));
    end
end
lambda4(kk,:) = eig(J4);
end
%% 3.5
t = 0:0.1:20;
n_0 = 0.25;
v_0 = -70;
I= @(t) 12;
y0 = [n_0,v_0];
[t,y] = ode15s(@(t,y)hhx2(t,y,I),t,y0);

plot(y(:,2),y(:,1),'color','magenta');
hold on;

n_0 = 0.45;
y0 = [n_0,v_0];
[t,y] = ode15s(@(t,y)hhx2(t,y,I),t,y0);

plot(y(:,2),y(:,1),'color','green');
legend("n-nullcline","v-nullcline","n0=0.25","n0=0.45");
hold off;
%%
function dy2dt = hhx(t,y,I)
dy2dt = zeros(4,1);
dy2dt(1) = ((2.5-0.1*(y(4)+74.44)+eps)/(exp(2.5-0.1*(y(4)+74.44))-1+eps))*(1-y(1))-(4*exp(-(y(4)+74.44)/18)+eps)*y(1); % dm/dt
dy2dt(2) = ((0.1-0.01*(y(4)+74.44)+eps)/(exp(1-0.1*(y(4)+74.44))-1+eps))*(1-y(2))-(0.125*exp(-(y(4)+74.44)/80)+eps)*y(2); % dn/dt
dy2dt(3) = (0.07*exp(-(y(4)+74.44)/20)+eps)*(1-y(3))-((1+eps)/(exp(3-0.1*(y(4)+74.44))+1+eps))*y(3); % dh/dt
dy2dt(4) = -(120*(y(1)^3)*y(3)*(y(4)-55)+36*(y(2)^4)*(y(4)+90)+0.3*(y(4)+60))+I(t); % dv/dt
end

function dy2dt = hhx2(t,y,I)
dy2dt = zeros(2,1);
dy2dt(1) = ((0.1-0.01*(y(2)+74.44)+eps)/(exp(1-0.1*(y(2)+74.44))-1+eps))*(1-y(1))-(0.125*exp(-(y(2)+74.44)/80)+eps)*y(1); % dn/dt
dy2dt(2) = -(120*((((2.5-0.1*(y(2)+74.44)+eps)/(exp(2.5-0.1*(y(2)+74.44))-1+eps))/((2.5-0.1*(y(2)+74.44)+eps)/(exp(2.5-0.1*(y(2)+74.44))-1+eps)+4*exp(-(y(2)+74.44)/18)+eps))^3)*(1-y(1))*(y(2)-55)+36*(y(1)^4)*(y(2)+90)+0.3*(y(2)+60))+I(t); % dv/dt
end



