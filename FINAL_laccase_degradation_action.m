ydata = ...
 [-0.12406 -0.58686 -0.68249];
xdata = ...
 [5 10 20];

fun = @(x, xdata) x(1) * xdata ./ (x(2) + xdata);
x0 = [-10 10];
x = lsqcurvefit(fun, x0, xdata, ydata);

times = linspace(xdata(1), xdata(end));
plot(xdata, ydata, '*r', times, fun(x,times), 'k', 'linewidth', 1.3)
legend('Measured data points','Fitted function')
xlabel('[Estrogen]_{initial} (µM)', 'fontsize', 14);
ylabel('Rate of reaction (\muM/min)', 'fontsize', 14);