close all
clear all

load('A2.mat')
titleSize = 16;
labelSize = 14;
legendSize = 12;
plotSize = 2;
markerSize = 8;
N = grids(1);
x = linspace(xl, xr, N)';
figure(3)
hold on
for i = 1:length(u_sol)
    plot(x, u_sol{i});
end
legend('eps=1','eps=0.1','eps=0.001','eps=0', 'FontSize', legendSize);
xlabel('x', 'FontSize', labelSize);
ylabel('u', 'FontSize', labelSize);
title('Solutions, different epsilon', 'FontSize', titleSize);
saveas(gcf, 'A2', 'png');