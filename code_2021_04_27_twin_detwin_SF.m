% This confirms that we can directly use the SF for detwinning 
addChenFunction;
eulerd = [15, 120, -150];   % assume euler angle, this is similar to the parent grain of Mg4Al_U2, ii=68, ID=75

stress = [-1 0 0; 0 0 0; 0 0 0];
m = angle2dcm(eulerd(1)/180*pi, eulerd(2)/180*pi, eulerd(3)/180*pi, 'zxz'); % transformation matrix

close all;

its = 3;    % twin system of interest

% [ref] look at twin system # 3
hcp_cell('euler', eulerd, 'stress',stress, 'ss', 24+its); 

ss = define_SS_cart('Mg','twin');

iss = 18 + its;

ng = ss(1,:,iss) * m;   % slip/twin plane normal in global coord
bg = ss(2,:,iss) * m;   % slip/twin direction in global coord
SF = ng * stress * bg';

eulerd_t = euler_by_twin(eulerd, its, 'Mg');   % twinned orientation

% [def] look at the same twin system #, in tension with -stress
hcp_cell('euler', eulerd_t, 'stress', -stress, 'ss', 24+its); 