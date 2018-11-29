I = importdata('./bigScan/initial_E_200.0_dr_0.06_mult_50._spec_0.out');
F = importdata('./bigScan/final_E_200.0_dr_0.06_mult_50._spec_0.out');

R = 1.67225;
halfMi = 1.6726219E-27 / 2 * 6.242E+18; %convert to eV

figure;
scatter(I(:, 2).^2 * halfMi, I(:, 1).^2 * halfMi, '.'); hold on;
% figure;
scatter(F(:, 2).^2 * halfMi, F(:, 1).^2 * halfMi,'.'); hold on;
plot(I(:, 2).^2 * halfMi, I(:, 2).^2.*(R-1) * halfMi);
title('200eV');
ylabel('E_{||} (eV)');
xlabel('E_{perp} (eV)');
legend('Initial', 'Lost');
set(gca, 'FontSize', 20);
% axis equal;
% axis equal;

% % lost = zeros(len(F200));
% % trapped = zeros(len(I200) - len(F200));

% % contourf(I200(:, 2).^2 * halfMi, I200(:, 1).^2 * halfMi, label);
% % figure;
% % heatmap(label, I200(:, 2).^2 * halfMi, I200(:, 1).^2 * halfMi)
% Eperp = I200(:, 2).^2 * halfMi;
% Epar = I200(:, 1).^2 * halfMi;
% [EPP, EPR] = meshgrid(Eperp, Epar);
% 
% onehot = ismember(I200, F200);
% label = diag(onehot(:, 1));
% figure;
% contourf(EPP, EPR, label)
