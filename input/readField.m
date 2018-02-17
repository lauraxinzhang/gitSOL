Field = importdata('LTX_Apr29_474-fields.dat',' ', 2);

R = reshape(Field.data(:,1), [260,260]);
Z = reshape(Field.data(:,2), [260,260]);

BR = reshape(Field.data(:,3), [260,260]);
BT = reshape(Field.data(:,4), [260,260]);
BZ = reshape(Field.data(:,5), [260,260]);

% figure;
% contour(R, Z, BT);

figure;
quiver(R, Z, BR, BZ)
streamline(R, Z, BR, BZ, [0.5, 0.51, 0.52, 0.53], [0, 0, 0, 0], [0.01, 100000])