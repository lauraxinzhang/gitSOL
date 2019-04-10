F = importdata('passingSolution.out', ',');

T = F(:, 1);
R = F(:, 2);
Phi = F(:, 3);

TM = reshape(T, [6, 20]);
RM = reshape(R, [6, 20]);
PhiM = reshape(Phi, [6, 20]);