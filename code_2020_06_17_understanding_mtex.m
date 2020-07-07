%% Code to understand Mtex 
addChenFunction;
cd('D:\p\m\mtex-5.3.1');
startup_mtex();

%% Methods to load EBSD data
clc
% path to files
pname = 'D:\WE43_T6_C1\EBSD Data';
% which files to be imported
fname = [pname '\WE43_T6_C1_Stitched_Dialated_rewrite.ang'];
cd(pname);

%% mtex way
% crystal symmetry. Note that the loaded .cif seems like to be Oxford convention. 
cs1 = crystalSymmetry('6/mmm', [3.2 3.2 5.2], 'X||a*', 'Y||b', 'Z||c', 'mineral', 'Mg');   % cs1 = crystalSymmetry.load('Mg-Magnesium.cif');
% cs2 is the one read from .ang file generated from Edax-OIM
cs2 = crystalSymmetry('622', [3.2 3.2 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Magnesium', 'color', [0.53 0.81 0.98]);
CS = {cs1, cs2};

% plotting convention
setMTEXpref('xAxisDirection','east');      % or 'east', ...
setMTEXpref('zAxisDirection','intoPlane'); % or 'intoPlane'

% Import the Data, create an EBSD variable containing the data. 
% The first data point should have euler angle about [11.5, 5.7, 2.9] degrees  
ebsd = EBSD.load(fname,CS,'interface','ang');   % can append flag (..., 'ang', 'convertEuler2SpatialReferenceFrame'/'convertEuler2SpatialReferenceFrame','setting 1/2/3/4')
euler_raw = [ebsd(1).rotations.phi1/pi*180, ebsd(1).rotations.Phi/pi*180, ebsd(1).rotations.phi2/pi*180]

% (0) Ground truth, using my method, for setting 2, align euler to sample
[a,b,c] = align_euler_to_sample(euler_raw(1), euler_raw(2), euler_raw(3), 0, -90, 180, 0);
euler_my_method = [a,b,c]

% (1) Can use flags to import --> This is the same as my treatment
ebsd = EBSD.load(fname,CS,'interface','ang','convertEuler2SpatialReferenceFrame','setting 2');    
euler_flag = [ebsd(1).rotations.phi1/pi*180, ebsd(1).rotations.Phi/pi*180, ebsd(1).rotations.phi2/pi*180]

% (2) A more flexable and explicit way: try to correct Data by (1) define a rot, and then (2) apply the rot
rot = rotation.byEuler(-90*degree,180*degree,0*degree); % This defines a rotation
ebsd = rotate(EBSD.load(fname,CS,'interface','ang'), rot, 'keepXY');       % (1) No flag rotate both, or (2) 'keepXY', (3) 'keepEuler', to rotate only one
euler_explicit = [ebsd(1).rotations.phi1/pi*180, ebsd(1).rotations.Phi/pi*180, ebsd(1).rotations.phi2/pi*180]


%% Ploting crystlas: (1) croping, indexing
% (1) My method, ground truth data
[phi1,phi,phi2,x,y,IQ,CI,Phase,Intensity,Fit] = read_ang(fname);
phi1 = phi1/pi*180;
phi = phi/pi*180;
phi2 = phi2/pi*180;

iR_min = 240;
iR_max = 280;
iC_min = 100;
iC_max = 150;
x = x(iR_min:iR_max, iC_min:iC_max);
y = y(iR_min:iR_max, iC_min:iC_max);
phi1 = phi1(iR_min:iR_max, iC_min:iC_max);
phi = phi(iR_min:iR_max, iC_min:iC_max);
phi2 = phi2(iR_min:iR_max, iC_min:iC_max);
[phi1, phi, phi2] = align_euler_to_sample(phi1, phi, phi2, 0, -90, 180, 0);

% for example, look at a selected point
pos = [x(15,15),y(15,15)];
euler = [phi1(15,15),phi(15,15),phi2(15,15)];
hcp_cell('euler',euler)
title(['euler = [',num2str(euler,'%8.1f'),']']);

%% (2) mtex, get orientation, and plot crystal
ebsd = EBSD.load(fname,CS,'interface','ang','convertEuler2SpatialReferenceFrame','setting 2');  
% cropping is by coordinate, not index
xmin = (iC_min-1) * 4;
xmax = (iC_max-1) * 4;
ymin = (iR_min-1) * 4;
ymax = (iR_max-1) * 4;
ebsd = ebsd(inpolygon(ebsd,[xmin, ymin, xmax-xmin, ymax-ymin]));

ori = ebsd.orientations;
plot(ebsd, ori)
% the same selected point, column major
pt_selected = ebsd(14*(iR_max-iR_min+1)+15)

% (*) Grain reconstruction
grains = calcGrains(ebsd, 'threshold', 10*degree)
% can smooth the grains to avoid the stair casing effect
% grains = smooth(grains,5);
hold on
plot(grains.boundary,'lineWidth',2)

% (*) Crystal shape
cS = crystalShape.hex(cs2)
% select only grains with more then 100 pixels
grains_selected = grains(grains.grainSize > 100);
plot(grains_selected, 0.7*cS, 'FaceColor',[0.3 0.5 0.3])

%% compute the ODF, and plot pole figure
odf = calcDensity(ori)
plotPDF(odf,Miller({0,0,0,1},{1,0,-1,0},cs2),'antipodal')
mtexColorbar;

%% Some concepts: rotation, orientation
close all;
clc;
% Edax convention, the euler angle defines how to rotation Lab-coordinate
% to align with Sample-coordinate by 3 succesive rotations about ZXZ
euler = [10, 5, 0]
% Plot in my default : x-right, y-down.
hcp_cell('euler',euler);

% Definition of rotation in MTEX
rot = rotation.byEuler(euler(1)*degree, euler(2)*degree, euler(3)*degree, 'Bunge')

% In MTEX, orientatin = rotation + crystal symmetry
cs = crystalSymmetry.load('Mg-Magnesium.cif');
ori = orientation.byEuler(euler(1)*degree, euler(2)*degree, euler(3)*degree, 'zxz', cs)

%% Quaternion expression, they are the SAME
% This is the built-in Matlab function, Q = [q0; q1,q2,q3], 
% = [cos(w/2); sin(w/2)*n], where (w,n) is the (angle, axis) expression
Q = angle2quat(euler(1)*degree, euler(2)*degree, euler(3)*degree,'zxz')

% The Mtex version of quaternion.
quaternion_mtex = quaternion(rot)

%% Transformation/Rotation matrix. 'Rotation matrix' is the SAME, where 'rotation matrix' is the transpose of 'transformation matrix'
% (1) Ground truth data, my method:
% Edax use 'passive rotation'
% m(i,j) = cos(x'_rotated, x_ref), so the row of the transformation matrix coincide with the new position of the rotated axis
transformation_matrix = angle2dcm(euler(1)/180*pi, euler(2)/180*pi, euler(3)/180*pi, 'zxz');
% the rotation matrix is the transpose
rotation_matrix = transformation_matrix'   

% (2) The mtex method is the same
% Another common way to represent rotations is by 3x3 matrices.
% The column of such a rotation matrix coincide with the new positions of the x, y and z vector after the rotation
rot.matrix % mtex convention, active rotation, rot.matrix = R = m'

%% uninstall mtex
cd('D:\p\m\mtex-5.3.1');
uninstall_mtex();
