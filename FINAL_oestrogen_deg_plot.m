y0 = 20;
p1 = -1.77;
p2 = 28.5;
f = @(t,x)(p1*x/(p2+x));
[t,y] = ode45(f, [1,100], y0);
plot(t, y, 'color', 'r', 'linewidth', 1.3);
xlabel('Time (minutes)', 'fontsize', 14);
ylabel('[Estrogen] (\muM)', 'fontsize', 14);