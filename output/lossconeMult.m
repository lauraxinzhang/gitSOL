I = importdata('./initial_Ti_200_Te_200_dr_0.06_mult_0.0_spec_0.out');
F = importdata('./final_Ti_200_Te_200_dr_0.06_mult_0.0_spec_0.out');

R = 1.67225;
Rminus = 1.72582;
Rplus = 1.62943;
halfMi = 1.6726219E-27 / 2 * 6.242E+18; %convert to eV

figure;
scatter(I(:, 2).^2 * halfMi, I(:, 1).^2 * halfMi, '.'); hold on;
% figure;
scatter(F(:, 2).^2 * halfMi, F(:, 1).^2 * halfMi,'.'); hold on;
plot(I(:, 2).^2 * halfMi, I(:, 2).^2.*(R-1) * halfMi);
plot(I(:, 2).^2 * halfMi, I(:, 2).^2.*(Rplus-1) * halfMi);
plot(I(:, 2).^2 * halfMi, I(:, 2).^2.*(Rminus-1) * halfMi);


title('200eV');
ylabel('E_{||} (eV)');
xlabel('E_{perp} (eV)');
legend('Initial', 'Lost');
set(gca, 'FontSize', 20);


I3 = importdata('./initial3_Ti_200_Te_200_dr_0.06_mult_0.0_spec_0.out');
F3 = importdata('./final3_Ti_200_Te_200_dr_0.06_mult_0.0_spec_0.out');

figure;
scatter3(I3(:,1), I3(:,2), I3(:, 3)); hold on;
scatter3(F3(:,1), F3(:,2), F3(:, 3), 'r');
figure;
scatter3(I3(:,1).^2 * halfMi, I3(:,2).^2 * halfMi, I3(:,3).^2 * halfMi); hold on;
scatter3(F3(:,1).^2 * halfMi, F3(:,2).^2 * halfMi, F3(:,3).^2 * halfMi, 'r');
xlabel('ER');
ylabel('Ephi');
zlabel('Ez');

