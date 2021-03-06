F = importdata('passingTest.out', ',');

R = F(:, 1);
Z = F(:, 2);
B = F(:, 3);

% RM = reshape(R, [259, 259]);
% ZM = reshape(Z, [259, 259]);
% BM = reshape(B, [259, 259]);

RM = reshape(R, [260, 260]);
ZM = reshape(Z, [260, 260]);
BM = reshape(B, [260, 260]);

figure;
set(gcf, 'renderer', 'zbuffer');
ax(1) = axes;
% surf(RM, ZM, BM, 'EdgeColor', 'none')
% contourf(RM, ZM, BM, 30);
mirror = contourf(RM, ZM, BM, 100,'edgecolor','none');
c1 = colorbar;
% c1.Ruler.Scale = 'log'
% caxis([0.5, 3.8745])
% set(c1, 'ylim', [-20 30])

% ylabel(c1, 'Electric Field (V/m)', 'FontSize', 16);
ylabel(c1, 'Normalized Potential (e\phi / Te)', 'FontSize', 16);

xlabel('Major Radius (m)', 'FontSize', 16);
ylabel('Cylindrical Height (m)', 'FontSize', 16);

% c1.label.string = 'Mirror Ratio';
% c.Label.String = 'Pastuknov Potential (e\phi / Te)';

% c.Label.String = 'Flux';

% title('Ti/Te = 0.2')
title('Electric Field')
set(gca, 'FontSize', 16);
axis equal;
colormap jet;
cmap = colormap;
% ncmap = sqrt(cmap);
% colormap(ncmap);

% ax(2) = axes;
% contour(RFLUX, ZFLUX, FFLUX, linspace(-0.0065, 0, 25), 'LineColor',[0.7 0.7 0.7],'LineStyle','--', 'LineWidth', 0.75);
% axis equal;
% 
% pos = get(ax(1), 'Position');
% set(ax(2), 'Visible', 'off');
% set(ax(2), 'Position', pos);
% 
% % ax(3) = axes;
% % contour(RFLUX, ZFLUX, FFLUX, linspace(-0.0042, -0.001, 17), 'LineColor',[0.5 0.5 0.5],'LineStyle','-', 'LineWidth', 0.75);
% % axis equal;
% % 
% % pos = get(ax(1), 'Position');
% % set(ax(3), 'Visible', 'off');
% % set(ax(3), 'Position', pos);
% 
% hold on;
% 
% FF = importdata('../input/Shell_Coordinates.csv', ',', 2);
% rShell = FF.data(:,1)* 0.0254; 
% zShell = FF.data(:,2)*0.0254;
% zminus = FF.data(:,4) * 0.0254;
% plot(rShell, zShell,'k', 'color',[191/255 191/255 191/255],'LineWidth', 2);
% plot(rShell, zminus,'k', 'color',[191/255 191/255 191/255],'LineWidth', 2);
