Ti_list = {'20','40','100'};
Te_list = {'100'};
mult_list = {'0.0','0.1','1.0'};

Ti = Ti_list{2};
Te = Te_list{1};
mult = mult_list{1};
spec = 0;

lsize = 3;     % line size
msize = 10;    % marker size
isize = 20;    % legend icon size


% Color array in RGB
% Note: poop yellow = [0.929, 0.694, 0.125]
%       poop red = [0.635, 0.078, 0.184]
%       washed blue = [0.3, 0.3, 1]

clr = [0.635, 0.078, 0.184;...
       0.929, 0.694, 0.125;...
       0, 0.7, 0];



range = 3* str2num(Ti) * (1 - spec) + 3 * str2num(Te) * spec;

subtitle = sprintf('\\fontsize{14}Ti = %s eV, Te = 100 eV, {\\phi} = %s{\\phi}_P', Ti, mult);
if length(Ti) == 2
    
    Ti = [Ti,'.'];
    
end



fileSTR = ['_Ti_',Ti,'_Te_100_dr_0.08_mult_',mult,'_spec_',num2str(spec),'.out'];
initial = ['initial',fileSTR];
final = ['final',fileSTR];

I = importdata(initial, ',',1);
II = I.data;
FI = importdata(final);

% R = 2.42048;
R = 2.13569;
Rl = 2.44586;
Rr = 1.98869;
halfM = (1.6726219E-27 * (1 - spec) + 9.1E-31 * spec ) / 2 * 6.242E+18; %convert to eV

figure(1)
hold on
plot(-1,-1,'.','Color',clr(1,:),'MarkerSize',isize)
plot(-1,-1,'.','Color',clr(2,:),'MarkerSize',isize)

% plot(II(:, 2).^2 * halfM, II(:, 1).^2 * halfM, '.', 'MarkerSize', msize, ...
%      'Color', clr(1,:));
plot(II(:, 1).^2 * halfM, II(:, 2).^2 * halfM,  '.', 'MarkerSize', msize, ...
     'Color', clr(1,:)); 
% figure;
% plot(FI(:, 2).^2 * halfM, FI(:, 1).^2 * halfM,'.', 'MarkerSize', msize, ...
%      'Color', clr(2,:));
plot(FI(:, 1).^2 * halfM, FI(:, 2).^2 * halfM, '.', 'MarkerSize', msize, ...
     'Color', clr(2,:));
plot([0, range], [0, range/(R-1)], 'LineWidth',lsize,...
     'LineStyle', '-','Color', clr(3,:));
plot([0, range], [0, range/(Rl-1)], 'LineWidth',lsize,...
     'LineStyle', '-','Color', clr(3,:));
plot([0, range], [0, range/(Rr-1)], 'LineWidth',lsize,...
     'LineStyle', '-','Color', clr(3,:));

xlabel('E_{||} (eV)');
ylabel('E_{\perp} (eV)');
legend('Trapped', 'Lost');

if spec == 1
    Title = 'Electron Loss Boundary';
else
    Title = 'Ion Loss Boundary';
end

set(gca, 'FontSize', 20);
title({Title; subtitle});
axis([0 range 0 range])
hold off
% axis equal;

