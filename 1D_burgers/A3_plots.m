close all
clear all

load('A3.mat')
titleSize = 16;
labelSize = 14;
legendSize = 12;
plotSize = 2;
markerSize = 8;
N = grids(1);
x = linspace(xl, xr, N)';
figure(3)
symbols = {'-x', '-o', '-d', '-s', '-*'};
hold on
for i = 1:length(u_sol)
    N = grids(i);
    x = linspace(xl, xr, N)';
    plot(x, u_sol{i}, symbols{i}, 'MarkerIndices', 1:10:length(u_sol{i}), 'MarkerSize', markerSize);
end
legend('N=41','N=81','N=161','N=321','N=641')
xlabel('x', 'FontSize', labelSize);
ylabel('u', 'FontSize', labelSize);
title('Solutions, epsilon grid dependent', 'FontSize', titleSize);
saveas(gcf, 'A3', 'png');