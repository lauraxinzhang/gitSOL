coord = importdata('coordRZ_E_100.000000_dr_0.070000spec0mult1.000000.out');

rCoord = coord(:,1);
phiCoord = coord(:,2);
zCoord = coord(:,3);

figure;
ax(1) = axes;
plot(rCoord, zCoord, '.');
xlabel('Major Radius (m)', 'FontSize', 16);
ylabel('Cylindrical Height (m)', 'FontSize', 16);
set(gca, 'FontSize', 16);
axis equal;

% coordXYZ = importdata('coordXYZ_E_200.000000_dr_0.083000spec0mult1.000000.out');
% xcoord = coordXYZ(:,1);
% ycoord = coordXYZ(:,2);
% zcoord = coordXYZ(:,3);
% 
% % figure;
% % plot(xcoord, ycoord, '.');
% % xlabel('Major Radius (m)', 'FontSize', 16);
% % ylabel('Cylindrical Height (m)', 'FontSize', 16);
% % 
% % set(gca, 'FontSize', 16);
% % axis equal;
% 
% figure;
% 
% scatter3(xcoord, ycoord, zcoord,'.');
% axis equal;

FLUX = importdata('flux.out');
RFLUX = reshape(FLUX(:, 1), [260, 260]);
ZFLUX = reshape(FLUX(:, 2), [260, 260]);
FFLUX = reshape(FLUX(:, 3), [260, 260]);

ax(2) = axes;
contour(RFLUX, ZFLUX, FFLUX, linspace(-0.0065, 0, 25), 'LineColor',[0.7 0.7 0.7],'LineStyle','--', 'LineWidth', 0.75);
axis equal;

pos = get(ax(1), 'Position');
set(ax(2), 'Visible', 'off');
set(ax(2), 'Position', pos);

limits_x = xlim(ax(2));
limits_y = ylim(ax(2));
xlim(ax(1),limits_x)
ylim(ax(1),limits_y)

hold on;

FF = importdata('../input/Shell_Coordinates.csv', ',', 2);
rShell = FF.data(:,1)* 0.0254; 
zShell = FF.data(:,2)*0.0254;
zminus = FF.data(:,4) * 0.0254;
plot(rShell, zShell,'k', 'color',[191/255 191/255 191/255],'LineWidth', 2);
plot(rShell, zminus,'k', 'color',[191/255 191/255 191/255],'LineWidth', 2);

