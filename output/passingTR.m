F = importdata('passingSolution.out', ' ');

T = F(:, 1);
R = F(:, 2);
Phi = F(:, 3);

TM = reshape(T, [20, 6]);
RM = reshape(R, [20, 6]);
PhiM = reshape(Phi, [20, 6]);

surf(TM, RM, PhiM)