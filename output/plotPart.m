A = importdata('12:17/coordRZ_E_40.000000_dr_0.070000spec0mult0.200000.out');
E = importdata('12:17/totalEnergy_E_40.000000_dr_0.070000spec0mult0.200000.out');
moment = importdata('12:17/Mu_E_40.000000_dr_0.070000spec0mult0.200000.out');

figure;
ax(1) = axes;
R = A(:,1);
Z = A(:,3);
plot(R, Z);
axis equal;

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


mu = moment(:,1);
vperp = moment(:, 2);
modB = moment(:, 3);
mu0 = mu(1); 

figure;
subplot(4, 1, 1);
plot(mu/mu0);
ylabel('magnetic moment')
subplot(4, 1, 2);
plot(E/E(1));
ylabel('energy')
subplot(4, 1, 3);
plot(R);
ylabel('R')
subplot(4, 1, 4);
plot(Z);
ylabel('Z')
xlabel('time steps')


% figure;
% subplot(4, 1, 1);
% plot(mu/mu0);
% ylabel(mu);
% subplot(4, 1, 2);
% plot(vperp.*vperp./(2*modB));
% subplot(4, 1, 3);
% plot(vperp);
% subplot(4, 1, 4);
% plot(modB);



