days = 0:7;
sp10 = [0 13088.802 12909.932 11916.074 13150.497 13235.365 17118.611 16584.274];
sp20 = [0 12983.868 13133.829 10544.882 15461.936 16429.047 17892.291 20510.374];

figure(1);
disp 'Figure 1'
scatter(days, sp10, 25,'b','*')
hold on
plot(days, sp20, 'r*')
h = lsline;
set(h(1), 'color', 'r', 'linewidth', 1.3);
set(h(2), 'color', 'b', 'linewidth', 1.3);
Lr = [ones(size(h(1).XData(:))), h(1).XData(:)] \ h(1).YData(:);
Slope_red = Lr(2)
Intercept_red = Lr(1)
Lb = [ones(size(h(2).XData(:))), h(2).XData(:)] \ h(2).YData(:);
Slope_blue = Lb(2)
Intercept_blue = Lb(1)
xlabel('Day', 'fontsize', 14);
ylabel('Protein Band Intensity (x 10^{4})', 'fontsize', 14);
set(gca,'YTick', 0:5000:25000)
set(gca,'YTickLabel',0:0.5:2.5)
legend ('Band Intensity SP 10', 'Band Intensity SP 20');

eq10 = [0 0.261028173400952 0.257460993656295 0.237640620610701 0.262258548278498 0.263951059099595 0.341394098595995 0.330737889487587];
eq20 = [0 0.258935489108864 0.261926140652938 0.210295127635713 0.308355257518788 0.32764221883167 0.356823492149113 0.409035560396619];

figure(2);
disp 'Figure 2'
scatter(days, eq10, 25, 'b', '*')
hold on
plot(days, eq20, 'r*')
g = lsline;
set(g(1), 'color', 'r', 'linewidth', 1.3);
set(g(2), 'color', 'b', 'linewidth', 1.3);
Lr = [ones(size(g(1).XData(:))), g(1).XData(:)] \ g(1).YData(:);
Slope_red = Lr(2)
Intercept_red = Lr(1)
Lb = [ones(size(g(2).XData(:))), g(2).XData(:)] \ g(2).YData(:);
Slope_blue = Lb(2)
Intercept_blue = Lb(1)
xlabel('Day', 'fontsize', 14);
ylabel('Amount of protein (\mug)', 'fontsize', 14);
set(gca, 'YTick', 0:0.1:0.5);
legend('Estimated quantity of SP 10', 'Estimated quantity of SP 20');