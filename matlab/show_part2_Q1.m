function [] = show_part2_Q1(t,y,I)
m = y(:,1);
n = y(:,2);
h = y(:,3);
V = y(:,4);

figure();
subplot(2,2,1)
plot(t, m)
title("m(t)")
xlabel("t[msec]")
ylabel("m")

subplot(2,2,2)
plot(t, n)
title("n(t)")
xlabel("t[msec]")
ylabel("n")

subplot(2,2,3)
plot(t, h)
title("h(t)")
xlabel("t[msec]")
ylabel("h")

subplot(2,2,4)
plot(t, V)
title("V(t)")
xlabel("t[msec]")
ylabel("V[mV]")

sgtitle_str = "$m, n, h$ and $V$ over time for an input of: $I = " + ...
              num2str(I) + " \cdot u(t)$";
sgtitle(sgtitle_str, "interpreter", "latex");
end

