Ti_list = {'20','40','100'};
Te_list = {'100'};
mult_list = {'0.0','0.1','0.5','1.0'};
iter_list = {'100000', '1000000','3000000'};
R_list = {3};


Ti = Ti_list{3};
Te = Te_list{1};
mult = mult_list{4};
iter = iter_list{2};
R = R_list{1};
spec = 1;

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

range = 6* str2num(Ti) * (1 - spec) + 6 * str2num(Te) * spec;

% subtitle = sprintf('\\fontsize{14}Ti = %s eV, Te = 100 eV, {\\phi} = %s{\\phi_P}, iter = %s', Ti, mult, iter);
subtitle = sprintf('\\fontsize{14}Ti = %s eV, Te = 100 eV, {\\phi} = %s{\\phi_P}', Ti, mult);
if length(Ti) == 2
    Ti = [Ti,'.'];    
end



fileSTR = ['_Ti_',Ti,'_Te_100_mult_',mult,'_R_',num2str(R),'_iter_',iter ,'_spec_',num2str(spec)];
initial = ['../output_new/initial',fileSTR,'.out'];
final = ['../output_new/final',fileSTR,'.out'];

II = importdata(initial);
FI = importdata(final);

halfM = (1.6726219E-27 * (1 - spec) + 9.1E-31 * spec ) / 2 * 6.242E+18; %convert to eV

figure()
hold on
plot(-1,-1,'.','Color',clr(1,:),'MarkerSize',isize)
plot(-1,-1,'.','Color',clr(2,:),'MarkerSize',isize)


plot(II(:, 1).^2 * halfM, II(:, 2).^2 * halfM,  '.', 'MarkerSize', msize, ...
     'Color', clr(1,:)); 
plot(FI(:, 1).^2 * halfM, FI(:, 2).^2 * halfM, '.', 'MarkerSize', msize, ...
     'Color', clr(2,:));
 
% plot expected loss regions 
plot([0, range], [0, range/(R-1)], 'LineWidth',lsize,...
     'LineStyle', '-.','Color', clr(3,:));
 
% shift = 456 * str2num(mult);
shift = 169 * str2num(mult);

if spec == 1
    plot([shift, range], [0, (range - shift)/(R - 1)],...
        'LineWidth',lsize,'LineStyle', '-','Color', clr(3,:));
else
    plot([0, range], [shift/(R-1), (range + shift)/(R - 1)],...
        'LineWidth',lsize,'LineStyle', '-','Color', clr(3,:));    
end

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
path = ['./bin/cone',fileSTR,'.eps'];
print(path, '-depsc', '-painters');
