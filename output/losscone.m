I10 = importdata('initial_E_10.00_dr_0.06.out');
F10 = importdata('final_E_10.00_dr_0.06.out');

I100 = importdata('initial_E_100.0_dr_0.06.out');
F100 = importdata('final_E_100.0_dr_0.06.out');




I500 = importdata('initial_E_1000._dr_0.06.out');
F500 = importdata('final_E_1000._dr_0.06.out');

I40 = importdata('initial_E_40.00_dr_0.06.out');
F40 = importdata('final_E_40.00_dr_0.06.out');

R = 1.67225;
halfMi = 1.6726219E-27 / 2 * 6.242E+18; %convert to eV

figure;
scatter(I10(:, 2).^2 * halfMi, I10(:, 1).^2 * halfMi, '.'); hold on;
% figure;
scatter(F10(:, 2).^2 * halfMi, F10(:,1).^2 * halfMi,'.'); hold on;
plot(I10(:, 2).^2 * halfMi, I10(:, 2).^2.*(R-1) * halfMi);
title('10eV');
ylabel('E_{||} (eV)');
xlabel('E_{perp} (eV)');
legend('Initial', 'Lost');
% axis equal;

figure;
scatter(I40(:, 2).^2 * halfMi, I40(:, 1).^2 * halfMi, '.'); hold on;
% figure;
scatter(F40(:, 2).^2 * halfMi, F40(:,1).^2 * halfMi,'.'); hold on;
plot(I40(:, 2).^2 * halfMi, I40(:, 2).^2.*(R-1) * halfMi);
title('40eV');
ylabel('E_{||} (eV)');
xlabel('E_{perp} (eV)');
legend('Initial', 'Lost');
% axis equal;

figure;
scatter(I100(:, 2).^2 * halfMi, I100(:, 1).^2 * halfMi,'.'); hold on; 
% figure;
scatter(F100(:, 2).^2 * halfMi, F100(:,1).^2 * halfMi,'.'); hold on;
plot(I100(:, 2).^2 * halfMi, I100(:, 2).^2.*(R-1) * halfMi);
title('100eV');
ylabel('E_{||} (eV)');
xlabel('E_{perp} (eV)');
legend('Initial', 'Lost');
% axis equal;

figure;
scatter(I500(:, 2).^2 * halfMi, I500(:, 1).^2 * halfMi,'.'); hold on; 
% figure;
scatter(F500(:, 2).^2 * halfMi, F500(:,1).^2 * halfMi,'.'); hold on;
plot(I500(:, 2).^2 * halfMi, I500(:, 2).^2.*(R-1) * halfMi);
title('1000eV');
ylabel('E_{||} (eV)');
xlabel('E_{perp} (eV)');
legend('Initial', 'Lost');
% axis equal;

% figure;
% histogram(I100(:,1)); hold on;
% histogram(F100(:, 1))
% % figure;
% histogram2(I10(:, 2).^2 , I10(:, 1).^2); hold on;
% % figure;
% histogram2(F10(:, 2).^2 , F10(:, 1).^2);
% 
% figure;
% histogram2(I100(:, 2).^2 , I100(:, 1).^2); hold on;
% % figure;
% histogram2(F100(:, 2).^2 , F100(:, 1).^2);
% 
% figure;
% histogram2(I500(:, 2).^2 , I500(:, 1).^2); hold on;
% % figure;
% histogram2(F500(:, 2).^2 , F500(:, 1).^2);


% figure;
% histogram(VV(:,1)); hold on; histogram(VV(:,2));
% % histogram(VV(:,1), VV(:,2));
% 
% figure;
% histogram(FF(:,1)); hold on; histogram(FF(:,2));
% % histogram(FF(:,1), FF(:,2))