function  plotVelocityInChannel(x,y,U,V,meshtype,nu,T)

quiver(x',y',U,V)     % velocity arrows
% imagesc([0,3],[0,1],U.^2 + V.^2)     % colormap
hold on
viscircles([0.5,0.5],0.1,'Color','k');
rectangle('Position',[0 0 3 1])

if nu > 0
    titleStr = sprintf('Velocity | mesh: %s | visc: %0.3f', meshtype, nu);
else
    titleStr = sprintf('Velocity | mesh: %s | visc: variable', meshtype);
end
title(titleStr, 'FontSize', 10);
xlim([0 3])
ylim([0 1])
daspect([1 1 1])

dim = [0.65 0.5 0.15 0.15];
str = sprintf('time: %0.2f', T);
annotation('textbox',dim,String=str,FitBoxToText='on', ...
    BackgroundColor='white', ...
    HorizontalAlignment='center', VerticalAlignment='middle');

fname = sprintf('plots/Channel/%s/nu_%0.3f/U/%f.png', meshtype, nu, T);
% saveas(gcf, fname, 'png','Resolution',300);
exportgraphics(gcf,fname,'Resolution',300)
hold off

end