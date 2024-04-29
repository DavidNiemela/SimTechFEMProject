
close all
clear all

load('A1.mat');
T = 0.4;
titleSize = 16;
labelSize = 14;
legendSize = 12;
plotSize = 2;
markerSize = 8;

c = 2;
u_exact = @(x,t, eps) c - tanh((x+0.5-c.*t)./(2*eps));
h = (xr-xl)./(grids-1);

for i = 1:length(err)-1
   q(i) = log(err(i)/err(i+1))/(grids(i)/grids(i+1));
end

figure(2)
loglog(h,err, '-o','LineWidth', plotSize);
hold on
% loglog(h, (err(1)/h(1)^2)*h.^2, '--','LineWidth', plotSize);
% loglog(h, (err(1)/h(1)^3)*h.^3, '--','LineWidth', plotSize);
loglog(h, (err(1)/h(1)^2.5)*h.^2.5, '--','LineWidth', plotSize);
hold off
legend('Errors', '2.5th order', 'FontSize', legendSize)
title('Convergence', 'FontSize', titleSize)
xlabel('h', 'FontSize', labelSize)
ylabel('errors', 'FontSize', labelSize)
saveas(gcf, 'A1_conv', 'png');

figure(3)
symbols = {'-x', '-o', '-d', '-s', '-*'};
hold on
for i = 1:length(u_sol)
    N = grids(i);
    x = linspace(xl, xr, N)';
    plot(x, u_sol{i}, symbols{i}, 'MarkerIndices', 1:10:length(u_sol{i}), 'MarkerSize', markerSize);
end
plot(x, u_exact(x,T,epsilons));
legend('N=41', 'N=81','N=161','N=321','N=641', 'Exact', 'FontSize', legendSize);
title('Solutions', 'FontSize', titleSize);
xlabel('x', 'FontSize', labelSize);
ylabel('u', 'FontSize', labelSize)
hold off
saveas(gcf, 'A1_sol', 'png');
