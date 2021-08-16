% this code checks the SF of T1 twins vs C1 twins
addChenFunction;
eulerd = [0, 10, 10];   % assume euler angle, this is similar to the parent grain of Mg4Al_U2, ii=68, ID=75

stress = [0 0 0; 0 0 0; 0 0 -1];    % z-axis compression
m = angle2dcm(eulerd(1)/180*pi, eulerd(2)/180*pi, eulerd(3)/180*pi, 'zxz'); % transformation matrix

ss = define_SS_cart('Ti','twin');


clear SF;
for its = 1:12
    iss = 24+its;
    ng = ss(1,:,iss) * m;   % slip/twin plane normal in global coord
    bg = ss(2,:,iss) * m;   % slip/twin direction in global coord
    SF(its,1) = ng * stress * bg';
end