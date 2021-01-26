%% EBSD data for in-situ test, UM134 Mg_C2, tested on 2021-01-15 (14 load steps, 2 cycles)
% Quick analysis for twin pct. Need to come back later.
% Tasks:
% (1) Add boundary to selected grains and generate new EBSD/grain data
% (2) Align Euler angles to sample reference frame.

%% setup
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-01-15 UM_134 Mg_C2 insitu EBSD';
save_dir = [working_dir, '\analysis'];
cd(working_dir);

%% reference, iE=0
iE = 0;
[EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C2 grain_file_type_1 iE=',num2str(iE),'.txt']));
[EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C2 grain_file_type_2 iE=',num2str(iE),'.txt']));

% find column
column_index_1 = find_variable_column_from_grain_file_header(EBSD_header_1, ...
    {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
column_index_2 = find_variable_column_from_grain_file_header(EBSD_header_2, ...
    {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});

% read type-2 grain file and get average info for grains
gID_0 = EBSD_data_2(:,column_index_2(1));
gPhi1_0 = EBSD_data_2(:,column_index_2(2));
gPhi_0 = EBSD_data_2(:,column_index_2(3));
gPhi2_0 = EBSD_data_2(:,column_index_2(4));
gCenterX_0 = EBSD_data_2(:,column_index_2(5));
gCenterY_0 = EBSD_data_2(:,column_index_2(6));
gEdge_0 = EBSD_data_2(:,column_index_2(10));

gNNeighbors = EBSD_data_2(:,column_index_2(7));
gDiameter = EBSD_data_2(:,column_index_2(8));
gArea = EBSD_data_2(:,column_index_2(9));
gEdge = EBSD_data_2(:,column_index_2(10));
gNeighbors = EBSD_data_2(:,(column_index_2(7)+1):(size(EBSD_data_2,2)));

% EBSD data, from type-1 grain file. (column, data) pair:
% (1,phi1) (2,phi) (3,phi2) (4,xMicron) (5,yMicron) (6,IQ) (7,CI) (8,Fit) (9,grain-Id) (10,edgeGrain?)
% Read EBSD data.  IQ,CI,Fit are not needed for now, but might need in future
x = EBSD_data_1(:,column_index_1(5));
y = EBSD_data_1(:,column_index_1(6));
unique_x = unique(x(:));
ebsdStepSize = unique_x(2) - unique_x(1);
mResize = (max(x(:)) - min(x(:)))/ebsdStepSize + 1;
nResize = (max(y(:)) - min(y(:)))/ebsdStepSize + 1;

x = reshape(EBSD_data_1(:,column_index_1(5)),mResize,nResize)';
y = reshape(EBSD_data_1(:,column_index_1(6)),mResize,nResize)';
ID_0 = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
boundary_0 = find_one_boundary_from_ID_matrix(ID_0);

%% Select deformed iE (parent grain file) to select control grains [No need to run every time]
for iE = 1:13
close all;

[EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C2 parent_grain_file_type_1 iE=',num2str(iE),'.txt']));
[EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C2 parent_grain_file_type_2 iE=',num2str(iE),'.txt']));

gID = EBSD_data_2(:,column_index_2(1));
gPhi1 = EBSD_data_2(:,column_index_2(2));
gPhi = EBSD_data_2(:,column_index_2(3));
gPhi2 = EBSD_data_2(:,column_index_2(4));
gCenterX = EBSD_data_2(:,column_index_2(5));
gCenterY = EBSD_data_2(:,column_index_2(6));
gEdge = EBSD_data_2(:,column_index_2(10));

gNNeighbors = EBSD_data_2(:,column_index_2(7));
gDiameter = EBSD_data_2(:,column_index_2(8));
gArea = EBSD_data_2(:,column_index_2(9));
gNeighbors = EBSD_data_2(:,(column_index_2(7)+1):(size(EBSD_data_2,2)));


ID = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
[boundary, boundaryID, neighborID, tripleTF, tripleID, indTriple, triIDs] = find_one_boundary_from_ID_matrix(ID);

phi1 = reshape(EBSD_data_1(:,column_index_1(2)),mResize,nResize)';
phi = reshape(EBSD_data_1(:,column_index_1(3)),mResize,nResize)';
phi2 = reshape(EBSD_data_1(:,column_index_1(4)),mResize,nResize)';
% change it to degrees, if necessary
if max(phi1(:))<7 && max(phi(:))<7 && max(phi2(:))<7
    phi1 = phi1*180/pi();
    phi = phi*180/pi();
    phi2 = phi2* 180/pi();
end

% plot to help selection
myplot(x, y, ID_0, grow_boundary(boundary_0)); 
title('ID, iE=0','fontweight','normal');
set(gca,'fontsize',18);
label_map_with_ID(x,y,ID_0,gcf, [61, 55, 388, 391],'r',24,10);  % illustrate selected control grains  
print(fullfile(save_dir, ['a_ID_map_iE=0.tif']),'-r300','-dtiff');

myplot(x, y, ID, grow_boundary(boundary)); 
title(['ID, iE=',num2str(iE)],'fontweight','normal');
set(gca,'fontsize',18);
print(fullfile(save_dir, ['a_ID_map_iE=',num2str(iE),'.tif']),'-r300','-dtiff');
% label_map_with_ID(x,y,ID,gcf,[18,26,205, 71],'r',24,10);    % illustrate selected control grains

%% Geotrans the iE=0 map to overlay on iE>0 map, to overlay the same grains [no need to run every time]
% (step-1) Select no less than 3 grain pairs for control points, and record
grain_pair{1} = [61, 59;
    55, 54;
    388, 384;
    391, 387];
grain_pair{2} = [61, 59;
    55, 58;
    388, 402;
    391, 406];
grain_pair{3} = [61, 57;
    55, 60;
    388, 406;
    391, 411];
grain_pair{4} = [61, 60;
    55, 61;
    388, 403;
    391, 407];
grain_pair{5} = [61, 59;
    55, 60;
    388, 393;
    391, 396];
grain_pair{6} = [61, 59;
    55, 57;
    388, 388;
    391, 391];
grain_pair{7} = [61, 59;
    55, 55;
    388, 381;
    391, 384];
grain_pair{8} = [61, 59;
    55, 58;
    388, 387;
    391, 391];
grain_pair{9} = [61, 57;
    55, 58;
    388, 390;
    391, 394];
grain_pair{10} = [61, 56;
    55, 60;
    388, 400;
    391, 405];
grain_pair{11} = [61, 59;
    55, 61;
    388, 393;
    391, 399];
grain_pair{12} = [61, 57;
    55, 58;
    388, 382;
    391, 386];
grain_pair{13} = [61, 60;
    55, 57;
    388, 378;
    391, 381];

% (step-2) Rough align using the selected control grains. The result is already decent 
g_0 = grain_pair{iE}(:,1);    % ref at iE=0
g_iE = grain_pair{iE}(:,2);    % iE > 0, considered as deformed
[~, loc_0] = ismember(g_0, gID_0);
[~, loc_iE] = ismember(g_iE, gID);
cpFrom = [gCenterX_0(loc_0),gCenterY_0(loc_0)];     % cpFrom is from ref image (iE=0)
cpTo = [gCenterX(loc_iE),gCenterY(loc_iE)];         % cpTo is from deformed image (iE>0)

% The tform is to transform the [ref @ iE=0] to [deformed iE>0]
tform = fitgeotrans(cpFrom, cpTo, 'affine');    % [x_0, y_0, 1] * tform.T = [x_iE, y_iE, 1]
ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');
boundary_0_to_iE = find_boundary_from_ID_matrix(ID_0_to_iE);

% myplot(x,y, ID_0_to_iE, boundary_0_to_iE); 
% set(gca,'fontsize',18);
% title(['(rough) affine transformed iE=0 -> iE=',num2str(iE)], 'fontweight','normal');



% (step-3) Fine align again, based on id_link and use all non-edge grains, and show result! 
[ID_new, id_link_additional, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);

% [[remove]] '0's from linked ids
ind_r = find(id_link(:,1)==0 | id_link(:,2)==0);
id_link(ind_r,:) = [];
inds = isnan(ID_0_to_iE)|(ID_new==0); % old ID is nan, or matched ID=0
ID_new(inds) = nan;

% [[remove]] grains that are shown in BOTH the col_2 of id_link and col_2 of in_link_additional  
ind = ismember(id_link(:,2), id_link_additional(:,2));
id_link(ind,:) = [];

% [[remove]] edge grains in id_link 
g_0 = id_link(:,1);
g_iE = id_link(:,2);    
[~, loc_0] = ismember(g_0, gID_0);
[~, loc_iE] = ismember(g_iE, gID);
ind = (gEdge_0(loc_0)==1) | (gEdge(loc_iE)==1);
g_0(ind) = [];
g_iE(ind) = [];

% use cps to transform
[~, loc_0] = ismember(g_0, gID_0);
[~, loc_iE] = ismember(g_iE, gID);
cpFrom = [gCenterX_0(loc_0),gCenterY_0(loc_0)]; % cpFrom is from deformed image (iE>0)
cpTo = [gCenterX(loc_iE),gCenterY(loc_iE)];   % cpTo is from ref image (iE=0)

% The tform is to transform the 'deformed' back into the 'ref'
tform = fitgeotrans(cpFrom, cpTo, 'affine');    % [x_ref, y_ref, 1] * tform.T = [x_def, y_def, 1]

ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');
boundary_0_to_iE = find_boundary_from_ID_matrix(ID_0_to_iE);

% myplot(x,y, ID_0_to_iE, boundary_0_to_iE); 
% set(gca,'fontsize',18);
% title(['(final) affine transformed iE=0 -> iE=',num2str(iE)], 'fontweight','normal');



% (step-4) try to match grains based on their spatial location, and show the matching result. 
% Transparent grains are the ones without a good match, so the parent grain need to be divided into more grains.     
[ID_new, id_link_additional, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);
% remove 0s
ind_r = find(id_link(:,1)==0 | id_link(:,2)==0);
id_link(ind_r,:) = [];
inds = isnan(ID_0_to_iE)|(ID_new==0); % old ID is nan, or matched ID=0
ID_new(inds) = nan;
% ---> just for illustration, I don't want to show the newly assigned ID
inds = ismember(ID_new, id_link_additional(:,1));
ID_new(inds) = nan;

myplot(x,y, ID_new, boundary_0_to_iE);
title(['iE=0 transformed with matched ID at iE=',num2str(iE)],'fontweight','normal');
set(gca,'fontsize',18);
print(fullfile(save_dir, ['a_tform_matched_ID_map_iE=',num2str(iE),'.tif']),'-r300','-dtiff');
% text(210,170,10000,'70<=>71','fontsize',18,'color','r')

end


%% align euler, and save .mat grain_file
close all;

% regular grain file
for iE = 0:13    
    
    [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C2 grain_file_type_1 iE=',num2str(iE),'.txt']));
    [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C2 grain_file_type_2 iE=',num2str(iE),'.txt']));
    % find column
    column_index_1 = find_variable_column_from_grain_file_header(EBSD_header_1, ...
        {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
    column_index_2 = find_variable_column_from_grain_file_header(EBSD_header_2, ...
        {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});
    
    gID = EBSD_data_2(:,column_index_2(1));
    gPhi1 = EBSD_data_2(:,column_index_2(2));
    gPhi = EBSD_data_2(:,column_index_2(3));
    gPhi2 = EBSD_data_2(:,column_index_2(4));
    gCenterX = EBSD_data_2(:,column_index_2(5));
    gCenterY = EBSD_data_2(:,column_index_2(6));
    gEdge = EBSD_data_2(:,column_index_2(10));
    
    gNNeighbors = EBSD_data_2(:,column_index_2(7));
    gDiameter = EBSD_data_2(:,column_index_2(8));
    gArea = EBSD_data_2(:,column_index_2(9));
    gNeighbors = EBSD_data_2(:,(column_index_2(7)+1):(size(EBSD_data_2,2)));
    
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
    
    [phi1, phi, phi2] = align_euler_to_sample(phi1, phi, phi2, 'non-mtex', 90,180,0);
    [gPhi1, gPhi, gPhi2] = align_euler_to_sample(gPhi1, gPhi, gPhi2, 'non-mtex', 90,180,0);
    
    euler_aligned_to_sample = 1;
    save(fullfile(save_dir, ['UM134_Mg_C2_grain_file_iE_',num2str(iE),'.mat']),'euler_aligned_to_sample','ID','phi1','phi','phi2','x','y',...
        'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');
    
end

% regular parent grain file
for iE = 0:13    
    
    [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C2 parent_grain_file_type_1 iE=',num2str(iE),'.txt']));
    [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C2 parent_grain_file_type_2 iE=',num2str(iE),'.txt']));
    % find column
    column_index_1 = find_variable_column_from_grain_file_header(EBSD_header_1, ...
        {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
    column_index_2 = find_variable_column_from_grain_file_header(EBSD_header_2, ...
        {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});
    
    gID = EBSD_data_2(:,column_index_2(1));
    gPhi1 = EBSD_data_2(:,column_index_2(2));
    gPhi = EBSD_data_2(:,column_index_2(3));
    gPhi2 = EBSD_data_2(:,column_index_2(4));
    gCenterX = EBSD_data_2(:,column_index_2(5));
    gCenterY = EBSD_data_2(:,column_index_2(6));
    gEdge = EBSD_data_2(:,column_index_2(10));
    
    gNNeighbors = EBSD_data_2(:,column_index_2(7));
    gDiameter = EBSD_data_2(:,column_index_2(8));
    gArea = EBSD_data_2(:,column_index_2(9));
    gNeighbors = EBSD_data_2(:,(column_index_2(7)+1):(size(EBSD_data_2,2)));
    
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
    
    [phi1, phi, phi2] = align_euler_to_sample(phi1, phi, phi2, 'non-mtex', 90,180,0);
    [gPhi1, gPhi, gPhi2] = align_euler_to_sample(gPhi1, gPhi, gPhi2, 'non-mtex', 90,180,0);
    
    euler_aligned_to_sample = 1;
    save(fullfile(save_dir, ['UM134_Mg_C2_modified_parent_grain_file_iE_',num2str(iE),'.mat']),'euler_aligned_to_sample','ID','phi1','phi','phi2','x','y',...
        'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');
    
end
