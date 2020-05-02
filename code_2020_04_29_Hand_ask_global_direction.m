
%% Keep track of the data
clc;
close all;
load('E:\Ti7Al_E1_insitu_tension\Analysis_by_Matlab\Ti7Al_E1_organized.mat', 'gPhi1','gPhi','gPhi2','gID');
ind = find(gID==226);
euler_226 = [gPhi1(ind), gPhi(ind), gPhi2(ind)];
disp('In organized data, euler angle for grain #2, ID: 226, is: ');
fprintf('[%.3f %.3f %.3f]\n', euler_226(1),euler_226(2),euler_226(3));

ind = find(gID==176);
euler_176 = [gPhi1(ind), gPhi(ind), gPhi2(ind)];
disp('In organized data, euler angle for grain #9, ID: 176, is: ');
fprintf('[%.3f %.3f %.3f]\n', euler_176(1),euler_176(2),euler_176(3));
%%

ss = define_SS('Mg','notwin');
ssc = define_SS_cart('Mg','notwin');

%% grain #2, ID=226, euler_organized = [-160.8850  162.4140  154.3300]
% ss #5: (10-10) [-12-10]  
% ss #4: (01-10) [2-1-10]
euler = euler_226;
m = angle2dcm(euler(1)/180*pi, euler(2)/180*pi, euler(3)/180*pi,'zxz');
disp(' ');
disp('grain #2, ID=226');

ssn = 5;
n0 = ss(1,:,ssn);
b0 = ss(2,:,ssn);
n =  ssc(1,:,ssn) * m;
b =  ssc(2,:,ssn) * m;
fprintf('\nss#%d: (%d %d %d %d)[%d %d %d %d]\n',ssn, n0(1),n0(2),n0(3),n0(4), b0(1),b0(2),b0(3),b0(4));
fprintf('In sample coordinate:\n slip plan normal direction = [%.3f %.3f %.3f]\n slip direction = [%.3f %.3f %.3f]\n', n(1),n(2),n(3),b(1),b(2),b(3));

ssn = 4;
n0 = ss(1,:,ssn);
b0 = ss(2,:,ssn);
n =  ssc(1,:,ssn) * m;
b =  ssc(2,:,ssn) * m;
fprintf('\nss#%d: (%d %d %d %d)[%d %d %d %d]\n',ssn, n0(1),n0(2),n0(3),n0(4), b0(1),b0(2),b0(3),b0(4));
fprintf('In sample coordinate:\n slip plan normal direction = [%.3f %.3f %.3f]\n slip direction = [%.3f %.3f %.3f]\n', n(1),n(2),n(3),b(1),b(2),b(3));

hcp_cell('euler',euler,'ss',[4,5])
%% grain #9, ID=176, euler_organized = [176.6880  132.5550  157.4800]
% ss #6:  (-1100) [-1-120]
% ss #11: (-1011) [-12-10]
% ss #4:  (01-10) [2-1-10]
% ss #1:  (0001) [2-1-10]

euler = euler_176;
disp(' ');
disp('grain #9, ID=176');

ssn = 6;
n0 = ss(1,:,ssn);
b0 = ss(2,:,ssn);
n =  ssc(1,:,ssn) * m;
b =  ssc(2,:,ssn) * m;
fprintf('\nss#%d: (%d %d %d %d)[%d %d %d %d]\n',ssn, n0(1),n0(2),n0(3),n0(4), b0(1),b0(2),b0(3),b0(4));
fprintf('In sample coordinate:\n slip plan normal direction = [%.3f %.3f %.3f]\n slip direction = [%.3f %.3f %.3f]\n', n(1),n(2),n(3),b(1),b(2),b(3));

ssn = 11;
n0 = ss(1,:,ssn);
b0 = ss(2,:,ssn);
n =  ssc(1,:,ssn) * m;
b =  ssc(2,:,ssn) * m;
fprintf('\nss#%d: (%d %d %d %d)[%d %d %d %d]\n',ssn, n0(1),n0(2),n0(3),n0(4), b0(1),b0(2),b0(3),b0(4));
fprintf('In sample coordinate:\n slip plan normal direction = [%.3f %.3f %.3f]\n slip direction = [%.3f %.3f %.3f]\n', n(1),n(2),n(3),b(1),b(2),b(3));

ssn = 4;
n0 = ss(1,:,ssn);
b0 = ss(2,:,ssn);
n =  ssc(1,:,ssn) * m;
b =  ssc(2,:,ssn) * m;
fprintf('\nss#%d: (%d %d %d %d)[%d %d %d %d]\n',ssn, n0(1),n0(2),n0(3),n0(4), b0(1),b0(2),b0(3),b0(4));
fprintf('In sample coordinate:\n slip plan normal direction = [%.3f %.3f %.3f]\n slip direction = [%.3f %.3f %.3f]\n', n(1),n(2),n(3),b(1),b(2),b(3));

ssn = 1;
n0 = ss(1,:,ssn);
b0 = ss(2,:,ssn);
n =  ssc(1,:,ssn) * m;
b =  ssc(2,:,ssn) * m;
fprintf('\nss#%d: (%d %d %d %d)[%d %d %d %d]\n',ssn, n0(1),n0(2),n0(3),n0(4), b0(1),b0(2),b0(3),b0(4));
fprintf('In sample coordinate:\n slip plan normal direction = [%.3f %.3f %.3f]\n slip direction = [%.3f %.3f %.3f]\n', n(1),n(2),n(3),b(1),b(2),b(3));

hcp_cell('euler',euler,'ss',[1,4,6,11])