
%% Make EBSD files for Reza:
addChenFunction;
p = 'E:\Mg4Al_3_sample_EBSD_Reza';
cd(p);

%% (1) Mg4Al_S1
f1 = 'Mg4Al_S1_grain_file_type_1.txt';
f2 = 'Mg4Al_S1_grain_file_type_2.txt';

% (1) load data by reading grain file
[EBSDdata1,EBSDheader1] = grain_file_read(fullfile(p,f1));
[EBSDdata2,EBSDheader2] = grain_file_read(fullfile(p,f2));
columnIndex1 = find_variable_column_from_grain_file_header(EBSDheader1,...
    {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
columnIndex2 = find_variable_column_from_grain_file_header(EBSDheader2,...
    {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});

% read type-2 grain file and get average info for grains
gID = EBSDdata2(:,columnIndex2(1));
gPhi1 = EBSDdata2(:,columnIndex2(2));
gPhi = EBSDdata2(:,columnIndex2(3));
gPhi2 = EBSDdata2(:,columnIndex2(4));

% EBSD data, from type-1 grain file. (column, data) pair:
% (1,phi1) (2,phi) (3,phi2) (4,xMicron) (5,yMicron) (6,IQ) (7,CI) (8,Fit) (9,grain-Id) (10,edgeGrain?)
% Read EBSD data.  IQ,CI,Fit are not needed for now, but might need in future
x = EBSDdata1(:,columnIndex1(5));
y = EBSDdata1(:,columnIndex1(6));
unique_x = unique(x(:));
ebsdStepSize = unique_x(2) - unique_x(1);
x_resize = (max(x(:)) - min(x(:)))/ebsdStepSize + 1;
y_resize = (max(y(:)) - min(y(:)))/ebsdStepSize + 1;

phi1 = reshape(EBSDdata1(:,columnIndex1(2)),x_resize,y_resize)';
phi = reshape(EBSDdata1(:,columnIndex1(3)),x_resize,y_resize)';
phi2 = reshape(EBSDdata1(:,columnIndex1(4)),x_resize,y_resize)';
% change it to degrees, if necessary
if max(phi1(:))<7 && max(phi(:))<7 && max(phi2(:))<7
    phi1 = phi1*180/pi();
    phi = phi*180/pi();
    phi2 = phi2* 180/pi();
end
x_um = reshape(EBSDdata1(:,columnIndex1(5)),x_resize,y_resize)';
y_um = reshape(EBSDdata1(:,columnIndex1(6)),x_resize,y_resize)';
ID = reshape(EBSDdata1(:,columnIndex1(1)),x_resize,y_resize)';

% crop, because forget to dewarp ...
if y_resize > x_resize
    disp(['x_dim = ',num2str(x_resize), ', y_dim = ',num2str(y_resize), ', so crop']);
    y_resize = x_resize;
    ID = ID(1:y_resize, 1:x_resize);
    x_um = x_um(1:y_resize, 1:x_resize);
    y_um = y_um(1:y_resize, 1:x_resize);
    phi1 = phi1(1:y_resize, 1:x_resize);
    phi = phi(1:y_resize, 1:x_resize);
    phi2 = phi2(1:y_resize, 1:x_resize);
end

% (2) Rotate euler angle
phiSys = [-90, 180, 0];
[phi1,phi,phi2] = align_euler_to_sample(phi1,phi,phi2,'none', phiSys(1),phiSys(2),phiSys(3)); % align euler angle to sample reference frame ------------ align.  UMich data is actually setting-1 !!!
% [q0,q1,q2,q3,phi1,phi,phi2] = regulate_euler_quat(phi1,phi,phi2);   % regulate the angles
[gPhi1,gPhi,gPhi2] = align_euler_to_sample(gPhi1,gPhi,gPhi2,'none', phiSys(1),phiSys(2),phiSys(3));

eulerAligned = 1;

save('Mg4Al_S1_EBSD_euler_aligned.mat','ID','x_um','y_um','gID','gPhi1','gPhi','gPhi2','phi1','phi','phi2','eulerAligned');

%% (2) Mg4Al_S2
f1 = 'Mg4Al_S2_grain_file_type_1.txt';
f2 = 'Mg4Al_S2_grain_file_type_2.txt';

% (1) load data by reading grain file
[EBSDdata1,EBSDheader1] = grain_file_read(fullfile(p,f1));
[EBSDdata2,EBSDheader2] = grain_file_read(fullfile(p,f2));
columnIndex1 = find_variable_column_from_grain_file_header(EBSDheader1,...
    {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
columnIndex2 = find_variable_column_from_grain_file_header(EBSDheader2,...
    {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});

% read type-2 grain file and get average info for grains
gID = EBSDdata2(:,columnIndex2(1));
gPhi1 = EBSDdata2(:,columnIndex2(2));
gPhi = EBSDdata2(:,columnIndex2(3));
gPhi2 = EBSDdata2(:,columnIndex2(4));

% EBSD data, from type-1 grain file. (column, data) pair:
% (1,phi1) (2,phi) (3,phi2) (4,xMicron) (5,yMicron) (6,IQ) (7,CI) (8,Fit) (9,grain-Id) (10,edgeGrain?)
% Read EBSD data.  IQ,CI,Fit are not needed for now, but might need in future
x = EBSDdata1(:,columnIndex1(5));
y = EBSDdata1(:,columnIndex1(6));
unique_x = unique(x(:));
ebsdStepSize = unique_x(2) - unique_x(1);
x_resize = (max(x(:)) - min(x(:)))/ebsdStepSize + 1;
y_resize = (max(y(:)) - min(y(:)))/ebsdStepSize + 1;

phi1 = reshape(EBSDdata1(:,columnIndex1(2)),x_resize,y_resize)';
phi = reshape(EBSDdata1(:,columnIndex1(3)),x_resize,y_resize)';
phi2 = reshape(EBSDdata1(:,columnIndex1(4)),x_resize,y_resize)';
% change it to degrees, if necessary
if max(phi1(:))<7 && max(phi(:))<7 && max(phi2(:))<7
    phi1 = phi1*180/pi();
    phi = phi*180/pi();
    phi2 = phi2* 180/pi();
end
x_um = reshape(EBSDdata1(:,columnIndex1(5)),x_resize,y_resize)';
y_um = reshape(EBSDdata1(:,columnIndex1(6)),x_resize,y_resize)';
ID = reshape(EBSDdata1(:,columnIndex1(1)),x_resize,y_resize)';

% crop, because forget to dewarp ...
if y_resize > x_resize
    disp(['x_dim = ',num2str(x_resize), ', y_dim = ',num2str(y_resize), ', so crop']);
    y_resize = x_resize;
    ID = ID(1:y_resize, 1:x_resize);
    x_um = x_um(1:y_resize, 1:x_resize);
    y_um = y_um(1:y_resize, 1:x_resize);
    phi1 = phi1(1:y_resize, 1:x_resize);
    phi = phi(1:y_resize, 1:x_resize);
    phi2 = phi2(1:y_resize, 1:x_resize);
end

% (2) Rotate euler angle
phiSys = [-90, 180, 0];
[phi1,phi,phi2] = align_euler_to_sample(phi1,phi,phi2,'none', phiSys(1),phiSys(2),phiSys(3)); % align euler angle to sample reference frame ------------ align.  UMich data is actually setting-1 !!!
% [q0,q1,q2,q3,phi1,phi,phi2] = regulate_euler_quat(phi1,phi,phi2);   % regulate the angles
[gPhi1,gPhi,gPhi2] = align_euler_to_sample(gPhi1,gPhi,gPhi2,'none', phiSys(1),phiSys(2),phiSys(3));

eulerAligned = 1;

save('Mg4Al_S2_EBSD_euler_aligned.mat','ID','x_um','y_um','gID','gPhi1','gPhi','gPhi2','phi1','phi','phi2','eulerAligned');

%% (3) Mg4Al_S3
f1 = 'Mg4Al_S3_grain_file_type_1.txt';
f2 = 'Mg4Al_S3_grain_file_type_2.txt';

% (1) load data by reading grain file
[EBSDdata1,EBSDheader1] = grain_file_read(fullfile(p,f1));
[EBSDdata2,EBSDheader2] = grain_file_read(fullfile(p,f2));
columnIndex1 = find_variable_column_from_grain_file_header(EBSDheader1,...
    {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
columnIndex2 = find_variable_column_from_grain_file_header(EBSDheader2,...
    {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});

% read type-2 grain file and get average info for grains
gID = EBSDdata2(:,columnIndex2(1));
gPhi1 = EBSDdata2(:,columnIndex2(2));
gPhi = EBSDdata2(:,columnIndex2(3));
gPhi2 = EBSDdata2(:,columnIndex2(4));

% EBSD data, from type-1 grain file. (column, data) pair:
% (1,phi1) (2,phi) (3,phi2) (4,xMicron) (5,yMicron) (6,IQ) (7,CI) (8,Fit) (9,grain-Id) (10,edgeGrain?)
% Read EBSD data.  IQ,CI,Fit are not needed for now, but might need in future
x = EBSDdata1(:,columnIndex1(5));
y = EBSDdata1(:,columnIndex1(6));
unique_x = unique(x(:));
ebsdStepSize = unique_x(2) - unique_x(1);
x_resize = (max(x(:)) - min(x(:)))/ebsdStepSize + 1;
y_resize = (max(y(:)) - min(y(:)))/ebsdStepSize + 1;

phi1 = reshape(EBSDdata1(:,columnIndex1(2)),x_resize,y_resize)';
phi = reshape(EBSDdata1(:,columnIndex1(3)),x_resize,y_resize)';
phi2 = reshape(EBSDdata1(:,columnIndex1(4)),x_resize,y_resize)';
% change it to degrees, if necessary
if max(phi1(:))<7 && max(phi(:))<7 && max(phi2(:))<7
    phi1 = phi1*180/pi();
    phi = phi*180/pi();
    phi2 = phi2* 180/pi();
end
x_um = reshape(EBSDdata1(:,columnIndex1(5)),x_resize,y_resize)';
y_um = reshape(EBSDdata1(:,columnIndex1(6)),x_resize,y_resize)';
ID = reshape(EBSDdata1(:,columnIndex1(1)),x_resize,y_resize)';

% crop, because forget to dewarp ...
if y_resize > x_resize
    disp(['x_dim = ',num2str(x_resize), ', y_dim = ',num2str(y_resize), ', so crop']);
    y_resize = x_resize;
    ID = ID(1:y_resize, 1:x_resize);
    x_um = x_um(1:y_resize, 1:x_resize);
    y_um = y_um(1:y_resize, 1:x_resize);
    phi1 = phi1(1:y_resize, 1:x_resize);
    phi = phi(1:y_resize, 1:x_resize);
    phi2 = phi2(1:y_resize, 1:x_resize);
end

% (2) Rotate euler angle
phiSys = [-90, 180, 0];
[phi1,phi,phi2] = align_euler_to_sample(phi1,phi,phi2,'none', phiSys(1),phiSys(2),phiSys(3)); % align euler angle to sample reference frame ------------ align.  UMich data is actually setting-1 !!!
% [q0,q1,q2,q3,phi1,phi,phi2] = regulate_euler_quat(phi1,phi,phi2);   % regulate the angles
[gPhi1,gPhi,gPhi2] = align_euler_to_sample(gPhi1,gPhi,gPhi2,'none', phiSys(1),phiSys(2),phiSys(3));

eulerAligned = 1;

save('Mg4Al_S3_EBSD_euler_aligned.mat','ID','x_um','y_um','gID','gPhi1','gPhi','gPhi2','phi1','phi','phi2','eulerAligned');

