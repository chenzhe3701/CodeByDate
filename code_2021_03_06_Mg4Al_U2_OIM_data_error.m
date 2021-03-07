 
% For sample Mg4Al_U2, the cleaned data for iE = 4 has something
% interesting. The average grain orientation of a grain was not corrected
% calculated by OIM analysis.

clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD';
save_dir = [working_dir, '\analysis'];
mkdir(save_dir);
sample_name = 'Mg4Al_U2';

%% Load data from grain file
iE = 4;
[EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,[sample_name,' grain_file_type_1 iE=',num2str(iE),'.txt']));
[EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,[sample_name,' grain_file_type_2 iE=',num2str(iE),'.txt']));
% find column
column_index_1 = find_variable_column_from_grain_file_header(EBSD_header_1, ...
    {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
column_index_2 = find_variable_column_from_grain_file_header(EBSD_header_2, ...
    {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});

gID = EBSD_data_2(:,column_index_2(1));
gPhi1 = EBSD_data_2(:,column_index_2(2));
gPhi = EBSD_data_2(:,column_index_2(3));
gPhi2 = EBSD_data_2(:,column_index_2(4));

x = EBSD_data_1(:,column_index_1(5));
y = EBSD_data_1(:,column_index_1(6));
unique_x = unique(x(:));
ebsdStepSize = unique_x(2) - unique_x(1);
mResize = (max(x(:)) - min(x(:)))/ebsdStepSize + 1;
nResize = (max(y(:)) - min(y(:)))/ebsdStepSize + 1;

ID = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
phi1 = reshape(EBSD_data_1(:,column_index_1(2)),mResize,nResize)';
phi = reshape(EBSD_data_1(:,column_index_1(3)),mResize,nResize)';
phi2 = reshape(EBSD_data_1(:,column_index_1(4)),mResize,nResize)';
x = reshape(EBSD_data_1(:,column_index_1(5)),mResize,nResize)';
y = reshape(EBSD_data_1(:,column_index_1(6)),mResize,nResize)';
% change it to degrees, if necessary
if max(phi1(:))<7 && max(phi(:))<7 && max(phi2(:))<7
    phi1 = phi1*180/pi();
    phi = phi*180/pi();
    phi2 = phi2* 180/pi();
end

indR_min = 491;
indR_max = 590;
indC_min = 1;
indC_max = 100;
ID_local = ID(indR_min:indR_max, indC_min:indC_max);
phi1_local = phi1(indR_min:indR_max, indC_min:indC_max);
phi_local = phi(indR_min:indR_max, indC_min:indC_max);
phi2_local = phi2(indR_min:indR_max, indC_min:indC_max);
x_local = x(indR_min:indR_max, indC_min:indC_max);
y_local = y(indR_min:indR_max, indC_min:indC_max);

% This can align euler
% [phi1, phi, phi2] = align_euler_to_sample(phi1, phi, phi2, 'non-mtex', 90,180,0);
% [gPhi1, gPhi, gPhi2] = align_euler_to_sample(gPhi1, gPhi, gPhi2, 'non-mtex', 90,180,0);

[nR,nC] = size(ID_local);
boundary_local = find_one_boundary_from_ID_matrix(ID_local);

%% generate IPF to illustrate
IPF_pixel = generate_IPF_map_pixel_wise(x_local, y_local, phi1_local, phi_local, phi2_local, boundary_local, [90,180,0], [1 0 0]);
IPF_grain = generate_IPF_map_grain_wise(ID_local, gID, gPhi1, gPhi, gPhi2, [90,180,0], [1 0 0], 1);

%% plot to illustrate error
close all;

% (1) label on pixel wise map
figure;
image(IPF_pixel);
ax = gca;
axis equal;
title('pixel wise IPF map in ED');

% point 1
indR = 20;
indC = 30;
euler_1 = [phi1_local(indR,indC), phi_local(indR,indC),phi2_local(indR,indC)];
hcp_cell('euler', euler_1, 'setting',1);
title('pt 1');
xy = [x_local(indR,indC),y_local(indR,indC)]
drawpoint(ax, 'position',[indC,indR], 'label',[' pt 1:, [', num2str(euler_1(1),4),' ,',num2str(euler_1(2),4),' ,',num2str(euler_1(3),4),']']);

% point 2
indR = 55;
indC = 15;
euler_2 = [phi1_local(indR,indC), phi_local(indR,indC),phi2_local(indR,indC)];
hcp_cell('euler', euler_2, 'setting',1);
title('pt 2');
xy = [x_local(indR,indC),y_local(indR,indC)]
drawpoint(ax, 'position',[indC,indR], 'label',[' pt 2:, [', num2str(euler_2(1),4),' ,',num2str(euler_2(2),4),' ,',num2str(euler_2(3),4),']']);

% (2) label on grain wise map
figure;
image(IPF_grain);
ax = gca;
axis equal;
title('grain wise IPF map in ED');

% point 1's grain 1
indR = 20;
indC = 30;
id_1 = ID_local(indR,indC); % 593
ind = find(gID==id_1);
euler_1 = [gPhi1(ind), gPhi(ind), gPhi2(ind)];
hcp_cell('euler', euler_1, 'setting',1);
title(['grain ',num2str(id_1)]);
drawpoint(ax, 'position',[indC,indR], 'label',[' g',num2str(id_1),':, [', num2str(euler_1(1),4),' ,',num2str(euler_1(2),4),' ,',num2str(euler_1(3),4),']']);

% point 2's grain 2
indR = 55;  % ylocal = 490+55-1=545
indC = 15;  % xlocal = 0+15-1=14
id_2 = ID_local(indR,indC);
ind = find(gID==id_2);
euler_1 = [gPhi1(ind), gPhi(ind), gPhi2(ind)];
hcp_cell('euler', euler_1, 'setting',1);
title(['grain ',num2str(id_2)]);
drawpoint(ax, 'position',[indC,indR], 'label',[' g',num2str(id_2),':, [', num2str(euler_1(1),4),' ,',num2str(euler_1(2),4),' ,',num2str(euler_1(3),4),']']);

%% now try to gather all data for the twinned grain 638, see if I can calculate the average orientation better 
dir_for_unit_cell = [save_dir,'\dir_for_unit_cell'];
mkdir(dir_for_unit_cell);

inds = (ID_local==638);
eulers = [phi1_local(inds),phi_local(inds),phi2_local(inds)];
xys = [x_local(inds),y_local(inds)];
ii = 1;
% hcp_cell('euler',eulers(1,:), 'setting', 1); title('ref orientation');
%% check ii=209 for pt 2
eulers_closest = zeros(size(eulers));
while ii <= size(eulers,1)
    % ===============> Find closest orientation among equivalent ones, and return euler
    eulers_closest(ii,:) = find_closest_orientation_hcp(eulers(ii,:), eulers(1,:));     
    
    % Plot 2 examples, showing two of the equivalent orientations
    if ismember(ii, [292, 293])
        hcp_cell('euler',eulers(ii,:), 'setting', 1); title('before correction');
        print(fullfile(dir_for_unit_cell, ['before_correction_',num2str(ii)]),'-dtiff');
        close;
        hcp_cell('euler', eulers_closest(ii,:), 'setting', 1); title('after correction');
        print(fullfile(dir_for_unit_cell, ['after_correction_',num2str(ii)]),'-dtiff');
        close;
    end
    
    ii = ii+1;
end
euler_avg = calculate_average_dominant_euler_hcp(eulers_closest)

%%














