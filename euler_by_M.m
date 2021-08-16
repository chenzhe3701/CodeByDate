% euler_d: euler angle (1x3) by degree
% M: transformation matrix
function euler_d_out = euler_by_M(euler_d_in, M_in)

% convert to radian
euler_r = euler_d_in/180*pi;

% convert euler_r to transformation matrix M
M = angle2dcm(euler_r(1),euler_r(2),euler_r(3),'zxz');

% multiply by input transformation
M = M_in * M;

% convert back to euler
[euler_r(1),euler_r(2),euler_r(3)] = dcm2angle(M,'zxz');

euler_d_out = euler_r/pi*180;

end