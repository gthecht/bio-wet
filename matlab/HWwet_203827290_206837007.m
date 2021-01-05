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
a_vec = linspace(a1, a2, 100);
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
%% 1.4

%% 2.1

%% 2.2

%% 2.3

%% 2.4

%% 2.5

%% 2.6

%% 2.7

%% 3.1

%% 3.2

%% 3.3

%% 3.4
