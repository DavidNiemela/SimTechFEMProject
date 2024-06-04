function plotPressureInLShape(p,t,P,meshtype,nu,T)

pdeplot(p,t,XYData=P,ColorMap="jet",Mesh="on")
hold on
pgon = polyshape([0.5,1,1,0,0,0.5],[0.5,0.5,0,0,1,1]);
plot(pgon, 'FaceColor','none')

if nu > 0
    titleStr = sprintf('Pressure | mesh: %s | visc: %0.3f', meshtype, nu);
else
    titleStr = sprintf('Pressure | mesh: %s | visc: variable', meshtype);
end
title(titleStr, 'FontSize', 10);
xlim([0 1])
ylim([0 1])
daspect([1 1 1])

dim = [0.62 0.7 0.15 0.15];
str = sprintf('time: %0.2f', T);
annotation('textbox',dim,String=str,FitBoxToText='on', ...
    BackgroundColor='white', ...
    HorizontalAlignment='center', VerticalAlignment='middle');

fname = sprintf('plots/LShape/%s/nu_%0.3f/P/%f.png', meshtype, nu, T);
% saveas(gcf, fname, 'png','Resolution',300);
exportgraphics(gcf,fname,'Resolution',300)
hold off

end