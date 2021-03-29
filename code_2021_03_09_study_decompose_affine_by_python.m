% Built rotation matrix, describe the rotation using passive rotation euler
% Example:
clc;
disp('Transformation: rotate +30 deg, then translate 2 in x');
euler_d = [30, 0, 0];
euler_r = euler_d/180*pi;
M = angle2dcm(euler_r(1), euler_r(2), euler_r(3), 'zxz');
R = M';  % we want R to perform active rotation
R2 = R(1:2, 1:2);

cpFrom = [0 0; 1 0; 0 1; 1 1]
cpTo = (R2 * cpFrom')'; % active rotation, [1 0] -> [0.86, 0.5]

% add a translation of 2 in x-directin to the transformation
cpTo(:,1) = cpTo(:,1) + 2

% Note that the description in Matlab is very very confusion.
% I fllowed this convention and it works well.
tform_by_fitgeotrans = fitgeotrans(cpFrom, cpTo, 'affine');
tform_by_maketform = maketform('affine', cpFrom(1:3,:), cpTo(1:3,:));

disp('tform.T = ')
disp(tform_by_fitgeotrans.T);

A2 = tform_by_fitgeotrans.T'

% For matlab, if use fitgeotrans(cpFrom(X,Y), cpTo(x,y), 'affine')
% then, tform convention is [X,Y,1] * tforms.T = [x,y,1]  
% Common and python transforms3d convention, is  A * [X;Y;1] = [x;y;1]

% decompose A into T,R,Z,S, using python function. 
% The input can be in different forms. They are equivalent.
[T,R,Z,S] = decompose_affine2d(A2)   % result = T,R,Z,S
[T,R,Z,S] = decompose_affine2d(tform_by_fitgeotrans) 
[T,R,Z,S] = decompose_affine2d(tform_by_maketform) 


%% Example, compose and decompose using python
T = py.numpy.array([20,30,40])  % translation
R = py.numpy.array([0 -1 0; 1 0 0; 0 0 1])  % rotation matrix
Z = py.numpy.array([2, 3, 4])   % zoom

A = py.transforms3d.affines.compose(T, R, Z)
result =py.transforms3d.affines.decompose(A)

TT = result(1)
RR = result(2)
ZZ = result(3)
SS = result(4)