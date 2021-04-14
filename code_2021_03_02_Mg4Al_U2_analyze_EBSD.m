
%% setup
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD';
save_dir = [working_dir, '\analysis'];
mkdir(save_dir);
% cd(working_dir);

sample_name = 'Mg4Al_U2';

%% Step-1: load .txt grain files, align_euler_to_sample and save as .mat

% (1.1) regular grain file
for iE = 0:13    
    
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
    save(fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']),'euler_aligned_to_sample','ID','phi1','phi','phi2','x','y',...
        'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');
    
end

% (1.2) parent grain file. [[[For iE=0, use grain_file as parent_grain_file]]]  
for iE = 1:13    
    [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,[sample_name,' parent_grain_file_type_1 iE=',num2str(iE),'.txt']));
    [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,[sample_name,' parent_grain_file_type_2 iE=',num2str(iE),'.txt']));
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
    save(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']),'euler_aligned_to_sample','ID','phi1','phi','phi2','x','y',...
        'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');
end
% For iE=0, use grain_file as parent_grain_file
copyfile(fullfile(save_dir, [sample_name,'_grain_file_iE_0.mat']), fullfile(save_dir, [sample_name,'_parent_grain_file_iE_0.mat']));

% save a copy after step 1
save_dir_1 = [save_dir,'\step-1'];
mkdir(save_dir_1);
for iE = 0:13
    copyfile(fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']), ...
           fullfile(save_dir_1, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));
    copyfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), ...
           fullfile(save_dir_1, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
end

%% Step-2: Analyze parent grain file: Select control points to affine transform and observe. Select grains to add grain boundaries. Update grain info for affected grains, and save data.
% Select deformed iE (parent grain file) to select control grains [No need to run every time]
for iE = 0:13
% iB=iE+1 for 0-based indexing, iB = iE for 1-based indexing
iB = iE+1;  
close all;    

% reference, iE=0
d = load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_0.mat']));

% read type-2 grain file and get average info for grains
gID_0 = d.gID;
gPhi1_0 = d.gPhi1;
gPhi_0 = d.gPhi;
gPhi2_0 = d.gPhi2;
gCenterX_0 = d.gCenterX;
gCenterY_0 = d.gCenterY;
gEdge_0 = d.gEdge;

gNNeighbors = d.gNNeighbors;
gDiameter = d.gDiameter;
gArea = d.gArea;
gNeighbors = d.gNeighbors;

x = d.x;
y = d.y;
ID_0 = d.ID;
boundary_0 = find_one_boundary_from_ID_matrix(ID_0);

% deformed iE
d = load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));

% read type-2 grain file and get average info for grains
gID = d.gID;
gPhi1 = d.gPhi1;
gPhi = d.gPhi;
gPhi2 = d.gPhi2;
gCenterX = d.gCenterX;
gCenterY = d.gCenterY;
gEdge = d.gEdge;
gNNeighbors = d.gNNeighbors;
gDiameter = d.gDiameter;
gArea = d.gArea;
gNeighbors = d.gNeighbors;

ID = d.ID;
phi1 = d.phi1;
phi = d.phi;
phi2 = d.phi2;
euler_aligned_to_sample = d.euler_aligned_to_sample;

[boundary, boundaryID, neighborID, tripleTF, tripleID, indTriple, triIDs] = find_one_boundary_from_ID_matrix(ID);

% plot to help selection
myplot(x, y, ID_0, grow_boundary(boundary_0)); 
title('ID, iE=0','fontweight','normal');
set(gca,'fontsize',18);
label_map_with_ID(x,y,ID_0,gcf, [18,31,108,120],'r',24,10);  % illustrate selected control grains  
print(fullfile(save_dir, ['a_ID_map_iE=0.tif']),'-r300','-dtiff');

myplot(x, y, ID, grow_boundary(boundary)); 
title(['ID, iE=',num2str(iE)],'fontweight','normal');
set(gca,'fontsize',18);
print(fullfile(save_dir, ['a_ID_map_iE=',num2str(iE),'.tif']),'-r300','-dtiff');
% label_map_with_ID(x,y,ID,gcf,[18,26,205, 71],'r',24,10);    % illustrate selected control grains

%% Geotrans the iE=0 map to overlay on iE>0 map, to overlay the same grains [no need to run every time]
% (step-1) Select no less than 3 grain pairs for control points, and record
grain_pair{1} = [18, 18;
    31, 31;
    108, 108;
    120, 120];
grain_pair{2} = [18, 17;
    31, 31;
    108, 105;
    120, 116];
grain_pair{3} = [18, 18;
    31, 29;
    108, 96;
    120, 110];
grain_pair{4} = [18, 19;
    31, 31;
    108, 98;
    120, 112];
grain_pair{5} = [18, 20;
    31, 33;
    108, 103;
    120, 117];
grain_pair{6} = [18, 19;
    31, 31;
    108, 95;
    120, 109];
grain_pair{7} = [18, 18;
    31, 33;
    108, 102;
    120, 116];
grain_pair{8} = [18, 18;
    31, 30;
    108, 99;
    120, 113];
grain_pair{9} = [18, 18;
    31, 32;
    108, 101;
    120, 114];
grain_pair{10} = [18, 18;
    31, 30;
    108, 99;
    120, 112];
grain_pair{11} = [18, 19;
    31, 30;
    108, 100;
    120, 114];
grain_pair{12} = [18, 19;
    31, 30;
    108, 101;
    120, 113];
grain_pair{13} = [18, 18;
    31, 31;
    108, 101;
    120, 114];
grain_pair{14} = [18, 18;
    31, 30;
    108, 100;
    120, 112];


% (step-2) Rough align using the selected control grains. The result is already decent 
g_0 = grain_pair{iB}(:,1);    % ref at iE=0
g_iE = grain_pair{iB}(:,2);    % iE > 0, considered as deformed

[~, loc_0] = ismember(g_0, gID_0);
[~, loc_iE] = ismember(g_iE, gID);
cpFrom = [gCenterX_0(loc_0),gCenterY_0(loc_0)];     % cpFrom is from ref image (iE=0)
cpTo = [gCenterX(loc_iE),gCenterY(loc_iE)];         % cpTo is from deformed image (iE>0)

% The tform is to transform the [ref @ iE=0] to [deformed iE>0]
tform = fitgeotrans(cpFrom, cpTo, 'affine');    % [x_0, y_0, 1] * tform.T = [x_iE, y_iE, 1]
ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');

% boundary_0_to_iE = find_boundary_from_ID_matrix(ID_0_to_iE);
% myplot(x,y, ID_0_to_iE, boundary_0_to_iE); 
% set(gca,'fontsize',18);
% title(['(rough) affine transformed iE=0 -> iE=',num2str(iE)], 'fontweight','normal');


% (step-3) Fine align again, based on id_link and use all non-edge grains, and show result! 
[ID_new, id_link_additional, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);

% [[remove]] '0's from linked ids
inds = find(id_link(:,1)==0 | id_link(:,2)==0);
id_link(inds,:) = [];
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

% (step-4) try to match grains based on their spatial location, and show the matching result. 
% Transparent grains are the ones without a good match, so the parent grain need to be divided into more grains.     
[ID_new, id_link_additional, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);
% remove 0s
ind = find(id_link(:,1)==0 | id_link(:,2)==0);
id_link(ind,:) = [];
inds = isnan(ID_0_to_iE)|(ID_new==0); % old ID is nan, or matched ID=0
ID_new(inds) = nan;
% ---> just for illustration, I don't want to show the newly assigned ID
inds = ismember(ID_new, id_link_additional(:,1));
ID_new(inds) = nan;

myplot(x,y, ID_new, boundary_0_to_iE);
title(['iE=0 transformed with matched ID at iE=',num2str(iE)],'fontweight','normal');
set(gca,'fontsize',18);
print(fullfile(save_dir, ['a_tform_matched_ID_map_iE=',num2str(iE),'.tif']),'-r300','-dtiff');



%% There might not be equal number of grains on the reference vs deformed maps. 
% (1) If there are more grains on the reference map than the deformed map,
% the extra grain on the reference map will not be able to find a match. If
% necessary, we want to select the grain on the deformed map, alter the
% tolerance angle to divide that grain into two (or more) grains. (2) If
% there are more grains on the deformed map, it's more difficult to check
% on the result 'deformed map with matched ID', but in this situation, the
% reference map might need to be processed, i.e., select a grain and reduce
% the tolerance angle to divide it into two grains.
% (3) In addition, sometimes we also need to merge grains on ref map

%% Record (1) the ID_list (on map_iE) to correct for each iE, and (2) the tolerance for each grain
% ID_list{iB=iE+1}. To treat 0-based indexing as 1 based indexing  
% ID_list{iB=iE+1} is the IDs need to be treated for iE = iB-1
ID_list{1} = [62,10]; 
tolerance_cell{1} = [3,3];
ID_list{2} = [2,16,50,87,127]; 
tolerance_cell{2} = [3,3,3,3,3];
ID_list{3} = [2,15,47,52,45, 25,80,117,116,121, 115,43,9]; 
tolerance_cell{3} = [3,3,3,3,3, 3,3,3,3,3, 3,3,3];
ID_list{4} = [16,49,53,85,27, 46,77,124,110]; 
tolerance_cell{4} = [3,3,3,3,3, 2,3,3,3];
ID_list{5} = [17,52,49,81,85, 115,125,129];
tolerance_cell{5} = [3,3,3,3,3, 3,3,3];
ID_list{6} = [10,16,47,50,54, 57,69,78,80,83, 121]; 
tolerance_cell{6} = [3,3,3,3,3, 3,3,2,3,3, 3];
ID_list{7} = [17,49,51,85,121, 126];
tolerance_cell{7} = [3,3,3,3,3, 3];
ID_list{8} = [2,9,17,47,52, 59,83,121,124]; 
tolerance_cell{8} = [3,3,3,3,2, 3,3,3,3];
ID_list{9} = [17,48,51,55,85, 122,125]; 
tolerance_cell{9} = [3,3,3,3,3, 3,3];
ID_list{10} = [2,16,47,50,54, 83,86,117,123]; 
tolerance_cell{10} = [3,3,3,3,3, 3,3,3,3];
ID_list{11} = [17,49,52,56,59, 81,84,112,125];
tolerance_cell{11} = [3,3,3,3,3, 3,3,3,3];
ID_list{12} = [10,16,47,52,54, 83,118,124]; 
tolerance_cell{12} = [3,3,3,3,3, 3,3,3];
ID_list{13} = [2,17,49,54,55, 83,118,119,124];
tolerance_cell{13} = [3,3,3,3,3, 3,3,3,3];
ID_list{14} = [10,17,46,51,83, 117,120,123]; 
tolerance_cell{14} = [3,3,3,3,3, 3,3,3];

% ID_merge_list{iB=iE+1} = [g1,->g2; ...]
ID_merge_list{1} = [];
ID_merge_list{2} = [];
ID_merge_list{3} = [];
ID_merge_list{4} = [84,85];
ID_merge_list{5} = [87,88];
ID_merge_list{6} = [];
ID_merge_list{7} = [];
ID_merge_list{8} = [];
ID_merge_list{9} = [];
ID_merge_list{10} = [];
ID_merge_list{11} = [];
ID_merge_list{12} = [];
ID_merge_list{13} = [];
ID_merge_list{14} = [];

% load the saved masks
try
    load(fullfile(save_dir,[sample_name,'_mask_cell.mat']),'mask_cell','ID_updated_cell');
catch
    mask_cell = cell(1,14);
end
 
iN = 1;
ID_updated = ID;
next_gID = max(gID) + 1;

for ii = 1:size(ID_merge_list{iB},1)
    ID_updated(ID_updated==ID_merge_list{iB}(ii,1)) = ID_merge_list{iB}(ii,2);
end
boundary = find_one_boundary_from_ID_matrix(ID_updated);   
boundary_new = boundary;
% In the next part, we need to check if the added boundary segment just
% create one addtional grain.

%% So, here, the function is to :
% (1) select the grain, (2) calculate local misorientation map, (3) find a
% line with the largest misorientation to divide the grain into two.

N = length(ID_list{iB});
while iN <= N
    try
        close(hf);
    end
    str = sprintf('iE=%d, iN=%d', iE, iN);
    disp(str);
    % example for debug, iE=1, grain 57
    ID_current = ID_list{iB}(iN);
    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
    indC_min = find(sum(ind_local, 1), 1, 'first');
    indC_max = find(sum(ind_local, 1), 1, 'last');
    indR_min = find(sum(ind_local, 2), 1, 'first');
    indR_max = find(sum(ind_local, 2), 1, 'last');
    
    ID_local = ID_updated(indR_min:indR_max, indC_min:indC_max);    % crop from ID_temp
    
    x_local = x(indR_min:indR_max, indC_min:indC_max);
    y_local = y(indR_min:indR_max, indC_min:indC_max);
    phi1_local = phi1(indR_min:indR_max, indC_min:indC_max);
    phi_local = phi(indR_min:indR_max, indC_min:indC_max);
    phi2_local = phi2(indR_min:indR_max, indC_min:indC_max);
    boundary_local = boundary_new(indR_min:indR_max, indC_min:indC_max);
    
    % for each pixel, find max misorientation within its four neighbors
    [nR,nC] = size(ID_local);
    misorientation_max = zeros(nR,nC);
    for iR = 1:nR
        for iC = 1:nC
            neighbor_orientation = [];
            % add upper neighbor, bottom neighbor, left neighbor, right neighbor
            if iR>1
                neighbor_orientation = [neighbor_orientation; phi1_local(iR-1,iC), phi_local(iR-1,iC), phi2_local(iR-1,iC)];
            end
            if iR<nR
                neighbor_orientation = [neighbor_orientation; phi1_local(iR+1,iC), phi_local(iR+1,iC), phi2_local(iR+1,iC)];
            end
            if iC>1
                neighbor_orientation = [neighbor_orientation; phi1_local(iR,iC-1), phi_local(iR,iC-1), phi2_local(iR,iC-1)];
            end
            if iC<nC
                neighbor_orientation = [neighbor_orientation; phi1_local(iR,iC+1), phi_local(iR,iC+1), phi2_local(iR,iC+1)];
            end
            
            euler_current = [phi1_local(iR,iC), phi_local(iR,iC), phi2_local(iR,iC)];
            
            misorientation = [];
            for ii = 1:size(neighbor_orientation,1)
                misorientation(ii) = calculate_misorientation_euler_d(euler_current, neighbor_orientation(ii,:), 'HCP');
            end
            misorientation_max(iR,iC) = max(misorientation);
        end
    end
    
    % if not previously processed, process, save the mask
    % plot the max misorientation map, use mask to select the region where you
    % want to add the misorientation boundary as a new boundary
    if (length(mask_cell{iB})>=iN) && ~isempty(mask_cell{iB}{iN})
        mask = mask_cell{iB}{iN};
    else
        hf2 = myplotm(boundary_local,'x',x_local,'y', y_local);
        label_map_with_ID(x_local,y_local,ID_local, gcf, ID_current, 'r');
        hf3= myplotm(misorientation_max);
        caxism([tolerance_cell{iB}(iN), 100]);
        h = drawpolygon;
        customWait(h);
        mask = h.createMask();
        mask_cell{iB}{iN} = mask;
    end
    
    tolerance = tolerance_cell{iB}(iN); % tolerance for this grain, default to 5
    boundary_local_new = (misorientation_max > tolerance)&(mask==1);
    boundary_local_new = double(ID_local~=ID_current | boundary_local_new);
    
    ID_local_new = find_ID_map_from_boundary_map(boundary_local_new);
    if length(unique(ID_local_new(:))) ~= 2
        % redraw
        mask_cell{iB}{iN} = [];     % delete mask
        myplot(boundary_local_new);
        warning(' More than one additional number of grains is created, you might want to do it again!');
    else
        % If split into g1 and g2, let g1 = 0, g2 = next_gID - ID_current, other grains = 0
        % ==> modification on 2021-02-25, need to check previous code for compatibility:  
        % Let the larger grain keep the ID_current.  
        if sum(ID_local_new(:)==1)/sum(ID_local(:)==ID_current) > 0.5
            ID_local_new(ID_local_new==1) = 0;
            ID_local_new(ID_local_new==2) = next_gID - ID_current;
        else
            ID_local_new(ID_local_new==2) = 0;
            ID_local_new(ID_local_new==1) = next_gID - ID_current;
        end
        ID_local_new(ID_local~=ID_current) = 0;
        ID_local_new = ID_local_new + ID_local;     % after adding, g2 will be next_gID
        ID_updated(indR_min:indR_max, indC_min:indC_max) = ID_local_new;    % update ID_new
        
        next_gID = next_gID + 1;    % update next_gID and iN
        iN = iN + 1;
        
        % update boundary, for illustration purpose.
        boundary_local_new = (misorientation_max > tolerance)&(mask==1);
        boundary_local_new = double(boundary_local | boundary_local_new);
        boundary_new(indR_min:indR_max, indC_min:indC_max) = boundary_local_new;
        try
            close(hf2);
            close(hf3);
        end
        hf = myplot(boundary_local_new);
    end
    
    ID_updated_cell{iB} = ID_updated;
    save(fullfile(save_dir, [sample_name,'_mask_cell.mat']),'mask_cell','ID_updated_cell');
end

close all;
myplot(boundary_new);

%% After adding more gb and get ID_updated, match ID_updated to ID_iE, update ID and grain info
% ID_updated: modified from the added grain boundaries
% ID: the target ID map at this iE

% umPerDp = 1;    % micron per data point in EBSD data
% data_in.umPerDp = umPerDp;
data_in.symmetry = 'hcp';
data_in.ID_target = ID;
data_in.ID_temp = ID_updated;
data_in.x = x;
data_in.y = y;
data_in.phi1 = phi1;
data_in.phi = phi;
data_in.phi2 = phi2;
data_in.gID = gID;
data_in.gPhi1 = gPhi1;
data_in.gPhi = gPhi;
data_in.gPhi2 = gPhi2;

data_out = match_ID_map_and_summarize(data_in);

ID = data_out.ID;   % the modified ID map
gID = data_out.gID;
gPhi1 = data_out.gPhi1;
gPhi = data_out.gPhi;
gPhi2 = data_out.gPhi2;
gCenterX = data_out.gCenterX;
gCenterY = data_out.gCenterY;
gNNeighbors = data_out.gNNeighbors;
gDiameter = data_out.gDiameter;
gArea = data_out.gArea;
gEdge = data_out.gEdge;
gNeighbors = data_out.gNeighbors;

save_dir_2 = [save_dir,'\step-2'];
mkdir(save_dir_2);

save(fullfile(save_dir_2, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']),'euler_aligned_to_sample','ID','phi1','phi','phi2','x','y',...
    'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');

close all;

disp(['finished iE=',num2str(iE)]);

% try to copy reference immediately. But if need to redo, need to copy from
% the saved file in step-1 folder.
if iE == 0
    copyfile(fullfile(save_dir_2, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), ...
        fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']) );
end

end


for iE = 0:13
    copyfile(fullfile(save_dir_2, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), ...
        fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']) );
end

%% special for this sample, modify 'grain file' at iE=5, divide twin grain into two grains.
for iE = 5
%% load 'grain file' data, plot map to select grain
iB = iE+1;

% deformed iE
d = load(fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));
% read type-2 grain file and get average info for grains
gID = d.gID;
gPhi1 = d.gPhi1;
gPhi = d.gPhi;
gPhi2 = d.gPhi2;

ID = d.ID;
phi1 = d.phi1;
phi = d.phi;
phi2 = d.phi2;
euler_aligned_to_sample = d.euler_aligned_to_sample;

boundary = find_one_boundary_from_ID_matrix(ID);

myplot(x, y, ID, grow_boundary(boundary)); 
title(['ID, iE=',num2str(iE)],'fontweight','normal');
set(gca,'fontsize',18);

%% record the grain to moidfy

ID_list{6} = [348]; 
tolerance_cell{6} = [3];

try
    load(fullfile(save_dir,[sample_name,'_twin_grain_mask_cell.mat']),'mask_cell','ID_updated_cell');
catch
    mask_cell = cell(1,14);
end

iN = 1;
ID_updated = ID;
next_gID = max(gID) + 1;

boundary = find_one_boundary_from_ID_matrix(ID_updated);   
boundary_new = boundary;

%% draw mask, select points as additional grain boundary to modify

N = length(ID_list{iB});
while iN <= N
    try
        close(hf);
    end
    str = sprintf('iE=%d, iN=%d', iE, iN);
    disp(str);
    % example for debug, iE=1, grain 57
    ID_current = ID_list{iB}(iN);
    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
    indC_min = find(sum(ind_local, 1), 1, 'first');
    indC_max = find(sum(ind_local, 1), 1, 'last');
    indR_min = find(sum(ind_local, 2), 1, 'first');
    indR_max = find(sum(ind_local, 2), 1, 'last');
    
    ID_local = ID_updated(indR_min:indR_max, indC_min:indC_max);    % crop from ID_temp
    
    x_local = x(indR_min:indR_max, indC_min:indC_max);
    y_local = y(indR_min:indR_max, indC_min:indC_max);
    phi1_local = phi1(indR_min:indR_max, indC_min:indC_max);
    phi_local = phi(indR_min:indR_max, indC_min:indC_max);
    phi2_local = phi2(indR_min:indR_max, indC_min:indC_max);
    boundary_local = boundary_new(indR_min:indR_max, indC_min:indC_max);
    
    % for each pixel, find max misorientation within its four neighbors
    [nR,nC] = size(ID_local);
    misorientation_max = zeros(nR,nC);
    for iR = 1:nR
        for iC = 1:nC
            neighbor_orientation = [];
            % add upper neighbor, bottom neighbor, left neighbor, right neighbor
            if iR>1
                neighbor_orientation = [neighbor_orientation; phi1_local(iR-1,iC), phi_local(iR-1,iC), phi2_local(iR-1,iC)];
            end
            if iR<nR
                neighbor_orientation = [neighbor_orientation; phi1_local(iR+1,iC), phi_local(iR+1,iC), phi2_local(iR+1,iC)];
            end
            if iC>1
                neighbor_orientation = [neighbor_orientation; phi1_local(iR,iC-1), phi_local(iR,iC-1), phi2_local(iR,iC-1)];
            end
            if iC<nC
                neighbor_orientation = [neighbor_orientation; phi1_local(iR,iC+1), phi_local(iR,iC+1), phi2_local(iR,iC+1)];
            end
            
            euler_current = [phi1_local(iR,iC), phi_local(iR,iC), phi2_local(iR,iC)];
            
            misorientation = [];
            for ii = 1:size(neighbor_orientation,1)
                misorientation(ii) = calculate_misorientation_euler_d(euler_current, neighbor_orientation(ii,:), 'HCP');
            end
            misorientation_max(iR,iC) = max(misorientation);
        end
    end
    
    % if not previously processed, process, save the mask
    % plot the max misorientation map, use mask to select the region where you
    % want to add the misorientation boundary as a new boundary
    if (length(mask_cell{iB})>=iN) && ~isempty(mask_cell{iB}{iN})
        mask = mask_cell{iB}{iN};
    else
        hf2 = myplotm(boundary_local,'x',x_local,'y', y_local);
        label_map_with_ID(x_local,y_local,ID_local, gcf, ID_current, 'r');
        hf3= myplotm(misorientation_max);
        caxism([tolerance_cell{iB}(iN), 100]);
        h = drawpolygon;
        customWait(h);
        mask = h.createMask();
        mask_cell{iB}{iN} = mask;
    end
    
    tolerance = tolerance_cell{iB}(iN); % tolerance for this grain, default to 5
    boundary_local_new = (misorientation_max > tolerance)&(mask==1);
    boundary_local_new = double(ID_local~=ID_current | boundary_local_new);
    
    ID_local_new = find_ID_map_from_boundary_map(boundary_local_new);
    if length(unique(ID_local_new(:))) ~= 2
        % redraw
        mask_cell{iB}{iN} = [];     % delete mask
        myplot(boundary_local_new);
        warning(' More than one additional number of grains is created, you might want to do it again!');
    else
        % If split into g1 and g2, let g1 = 0, g2 = next_gID - ID_current, other grains = 0
        % ==> modification on 2021-02-25, need to check previous code for compatibility:  
        % Let the larger grain keep the ID_current.  
        if sum(ID_local_new(:)==1)/sum(ID_local(:)==ID_current) > 0.5
            ID_local_new(ID_local_new==1) = 0;
            ID_local_new(ID_local_new==2) = next_gID - ID_current;
        else
            ID_local_new(ID_local_new==2) = 0;
            ID_local_new(ID_local_new==1) = next_gID - ID_current;
        end
        ID_local_new(ID_local~=ID_current) = 0;
        ID_local_new = ID_local_new + ID_local;     % after adding, g2 will be next_gID
        ID_updated(indR_min:indR_max, indC_min:indC_max) = ID_local_new;    % update ID_new
        
        next_gID = next_gID + 1;    % update next_gID and iN
        iN = iN + 1;
        
        % update boundary, for illustration purpose.
        boundary_local_new = (misorientation_max > tolerance)&(mask==1);
        boundary_local_new = double(boundary_local | boundary_local_new);
        boundary_new(indR_min:indR_max, indC_min:indC_max) = boundary_local_new;
        try
            close(hf2);
            close(hf3);
        end
        hf = myplot(boundary_local_new);
    end
    
    ID_updated_cell{iB} = ID_updated;
    save(fullfile(save_dir, [sample_name,'_twin_grain_mask_cell.mat']),'mask_cell','ID_updated_cell');
end

close all;
myplot(boundary_new);

%% After adding more gb and get ID_updated, match ID_updated to ID_iE, update ID and grain info
% ID_updated: modified from the added grain boundaries
% ID: the target ID map at this iE

% umPerDp = 1;    % micron per data point in EBSD data
% data_in.umPerDp = umPerDp;
data_in.symmetry = 'hcp';
data_in.ID_target = ID;
data_in.ID_temp = ID_updated;
data_in.x = x;
data_in.y = y;
data_in.phi1 = phi1;
data_in.phi = phi;
data_in.phi2 = phi2;
data_in.gID = gID;
data_in.gPhi1 = gPhi1;
data_in.gPhi = gPhi;
data_in.gPhi2 = gPhi2;

data_out = match_ID_map_and_summarize(data_in);

ID = data_out.ID;   % the modified ID map
gID = data_out.gID;
gPhi1 = data_out.gPhi1;
gPhi = data_out.gPhi;
gPhi2 = data_out.gPhi2;
gCenterX = data_out.gCenterX;
gCenterY = data_out.gCenterY;
gNNeighbors = data_out.gNNeighbors;
gDiameter = data_out.gDiameter;
gArea = data_out.gArea;
gEdge = data_out.gEdge;
gNeighbors = data_out.gNeighbors;

save_dir_2 = [save_dir,'\step-2'];
mkdir(save_dir_2);

save(fullfile(save_dir_2, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']),'euler_aligned_to_sample','ID','phi1','phi','phi2','x','y',...
    'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');

close all;

disp(['finished twin grain file at iE=',num2str(iE)]);

end

%% copy to analysis folder
for iE = 0:13
    try
        copyfile(fullfile(save_dir_2, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']), ...
            fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']) );
    end
end

%% Part-3: Affine transform ID map. Link grains, and modify linked grains to the same ID#.
% Generate geotrans/tform information at iE = 1 to 13, without showing results
% Link ids (save in 'tbl') at differnt iEs after loading the saved geotrans/tform. 
%   So, later, we can have grain files output again with different IDs.

% reference, iE=0
d = load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_0.mat']));

x = d.x;
y = d.y;
ID_0 = d.ID;
boundary_0 = find_one_boundary_from_ID_matrix(ID_0);
gID_0 = d.gID;
gPhi1_0 = d.gPhi1;
gPhi_0 = d.gPhi;
gPhi2_0 = d.gPhi2;
gCenterX_0 = d.gCenterX;
gCenterY_0 = d.gCenterY;
gEdge_0 = d.gEdge;

for iE = 1:13
    iB = iE + 1;
    disp(['iE=',num2str(iE)]);
    % load the modified EBSD data at iE
    load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), 'ID','gID','gEdge','gCenterX','gCenterY');
    
    % match IDs
    g_0 = grain_pair{iB}(:,1);    % ref at iE=0
    g_iE = grain_pair{iB}(:,2);    % iE > 0, considered as deformed
    [~, loc_0] = ismember(g_0, gID_0);
    [~, loc_iE] = ismember(g_iE, gID);
    cpFrom = [gCenterX_0(loc_0),gCenterY_0(loc_0)];     % cpFrom is from ref image (iE=0)
    cpTo = [gCenterX(loc_iE),gCenterY(loc_iE)];         % cpTo is from deformed image (iE>0)
    
    % [1st] rough tform. The tform is to transform the [ref @ iE=0] to [deformed iE>0]
    tform = fitgeotrans(cpFrom, cpTo, 'affine');    % [x_0, y_0, 1] * tform.T = [x_iE, y_iE, 1]
    tforms_0{iB} = tform;     % save the tform
    
    ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');
    boundary_0_to_iE = find_boundary_from_ID_matrix(ID_0_to_iE);
    
    [ID_new, id_link_additional, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);
    
    % [[remove]] '0's from linked ids
    ind_r = find(id_link(:,1)==0 | id_link(:,2)==0);
    id_link(ind_r,:) = [];
    inds = isnan(ID_0_to_iE)|(ID_new==0); % old ID is nan, or matched ID=0
    ID_new(inds) = nan;
    
    % Plot rough tformed map
    myplot(x,y, ID_new, boundary_0_to_iE);
    title(['iE=0 rough transformed with matched ID at iE=',num2str(iE)],'fontweight','normal');
    set(gca,'fontsize',18);
    print(fullfile(save_dir, ['b_rough_tform_matched_ID_map_iE=',num2str(iE),'.tif']),'-r300','-dtiff');
    close;
    
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
    
    % [2nd] fine tform.
    tform = fitgeotrans(cpFrom, cpTo, 'affine');    % [x_ref, y_ref, 1] * tform.T = [x_def, y_def, 1]
    tforms{iB} = tform;     % save the tform
    
    ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');
    boundary_0_to_iE = find_boundary_from_ID_matrix(ID_0_to_iE);
    
    [ID_new, id_link_additional, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);
    % remove 0s
    ind_r = find(id_link(:,1)==0 | id_link(:,2)==0);
    id_link(ind_r,:) = [];
    inds = isnan(ID_0_to_iE)|(ID_new==0); % old ID is nan, or matched ID=0
    ID_new(inds) = nan;
    
    % plot fine tformed map
    myplot(x,y, ID_new, boundary_0_to_iE);
    title(['iE=0 transformed with matched ID at iE=',num2str(iE)],'fontweight','normal');
    set(gca,'fontsize',18);
    print(fullfile(save_dir, ['b_tform_matched_ID_map_iE=',num2str(iE),'.tif']),'-r300','-dtiff');
    close;
    
    tbl_iE = array2table(id_link);
    tbl_iE.Properties.VariableNames = {'iE_0',['iE_',num2str(iE)]};
    
    if iE==1
        tbl = tbl_iE;
    else
        tbl = innerjoin(tbl,tbl_iE);
    end
end

save(fullfile(save_dir,'geotrans_and_id_link.mat'),'tforms_0','tforms','tbl');

%% part-4, modify linked ID numbers to be the same throughout iEs
load(fullfile(save_dir,'geotrans_and_id_link.mat'),'tforms','tbl');

maxID = 0;
for iE = 0:13
    f = matfile(fullfile(save_dir,[sample_name,'_parent_grain_file_iE_',num2str(iE)]));

    gID = f.gID;
    maxID = max(maxID, max(gID));
end
id_to_add = 10^(ceil(log10(maxID)));    % round to 1000 etc for adjustment. 

for iE = 0:13
    clear euler_aligned_to_sample;
    load(fullfile(save_dir,[sample_name,'_parent_grain_file_iE_',num2str(iE)]));

    if ~exist('euler_aligned_to_sample','var')
        warning('suggest to align_euler_to_sample first');
        euler_aligned_to_sample = 0;
    end
    
    % adjust, add a big enough value to prevent old unmodified grains having the same id as new modified grains  
    ID(ID>0) = ID(ID>0) + id_to_add;
    gID(gID>0) = gID(gID>0) + id_to_add;
    gNeighbors(gNeighbors>0) = gNeighbors(gNeighbors>0) + id_to_add;    % add
    
    id_link = tbl.Variables;
    g_list_target = id_link(:,1);
    g_list_source = id_link(:,iE+1) + id_to_add;  % Temporarily change the gid in the source gID list.  Remember to use 'iE+1'
    for ii = 1:length(g_list_target)
        id_source = g_list_source(ii);
        id_target = g_list_target(ii);
        % modify
        ID(ID==id_source) = id_target;
        
        ind = find(gID_0==id_target);
        eulerd_ref = [gPhi1_0(ind), gPhi_0(ind), gPhi2_0(ind)];
        ind = find(gID==id_source);
        eulerd = [gPhi1(ind), gPhi(ind), gPhi2(ind)];
        
        % eulerd = find_closest_orientation_hcp(eulerd, eulerd_ref);      % ==> Find Closest Orientation
        
        gPhi1(gID==id_source) = eulerd(1);
        gPhi(gID==id_source) = eulerd(2);
        gPhi2(gID==id_source) = eulerd(3);
        
        gID(gID==id_source) = id_target;
        gNeighbors(gNeighbors==id_source) = id_target;
    end
    
    boundary = find_boundary_from_ID_matrix(ID);
    myplot(x,y,mod(ID,10),boundary);
    clim = caxis;
%     caxis([0,maxID]);
    title(['iE=',num2str(iE),', mod(ID spatially aligned to iE=0,10) #']);
    print(fullfile(save_dir, ['matched_ID_map_iE=',num2str(iE),'.tif']),'-r300','-dtiff');
    close;
    
    save(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), 'euler_aligned_to_sample', ...
        'ID','phi1','phi','phi2','x','y',...
        'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');

end

save_dir_4 = [save_dir,'\step-4'];
mkdir(save_dir_4);
for iE = 0:13
    copyfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), ...
           fullfile(save_dir_4, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
end

%% part-5: find out variants
po_tolerance_angle = 10; % if child grain has misorientation < po_tolerance_angle with undeformed parent grain, it is considered as having parent orientation
twin_tolerance_angle = 10;  % if child grain has misorientation < twin_tolerance_angle to a potential twin variant, it is identified as that twin variant  
for iE = 1:13
    d = load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_0.mat']));
    gID_0 = d.gID;
    gPhi1_0 = d.gPhi1;
    gPhi_0 = d.gPhi;
    gPhi2_0 = d.gPhi2;
    ID_0 = d.ID;
    boundary_0 = find_boundary_from_ID_matrix(ID_0);
    
    d = load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    x = d.x;
    y = d.y;
    gID_p = d.gID;
    gPhi1_p = d.gPhi1;
    gPhi_p = d.gPhi;
    gPhi2_p = d.gPhi2;
    ID_p = d.ID;
    boundary_p = find_boundary_from_ID_matrix(ID_p);
        
    % data where twins (children) are individially labeled with IDs
    d = load(fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));
    gID_c = d.gID;
    gPhi1_c = d.gPhi1;
    gPhi_c = d.gPhi;
    gPhi2_c = d.gPhi2;
    ID_c = d.ID;
    phi1_c = d.phi1;
    phi_c = d.phi;
    phi2_c = d.phi2;
    
    ID_variant_grain_wise = zeros(size(ID_p));
    ID_variant_point_wise = zeros(size(ID_p));
    Misorientation_point_wise = zeros(size(ID_p));

    for ii = 1: length(gID_p)        
        id_p = gID_p(ii);
        disp(['current ID = ',num2str(id_p)]);
        
        % Find the parent orientation at iE = 0    
        ind_0 = find(gID_0 == id_p);
        euler_0 = [gPhi1_0(ind_0), gPhi_0(ind_0), gPhi2_0(ind_0)];
        
        if ~isempty(ind_0)
            % find out all the id numbers on ID_c overlap with the parent grain on ID_p
            id_c = unique(ID_c(ID_p == id_p));
            
            % Additionally, the child should have > 70% of its pixels overlap with the potential parent grain
            area_c = zeros(size(id_c));
            overlap_c = zeros(size(id_c));
            for kk = 1:length(area_c)
                area_c(kk) = sum(ID_c(:)==id_c(kk));
                overlap_c(kk) = sum((ID_p(:)==id_p)&(ID_c(:)==id_c(kk)));
            end
            overlap_pct_c = overlap_c./area_c;
            id_c(overlap_pct_c<0.75) = [];
            
            % Modify all the child orientation to crystallographically equivalent one that is closest to euler_0  
            misorientation = [];
            col_ID = [];
            col_euler = [];
            for jj = 1:length(id_c)
                id = id_c(jj);
                ind = (gID_c == id);
                euler_id = [gPhi1_c(ind), gPhi_c(ind), gPhi2_c(ind)];
                
                % =================> Important, here modified the grain info of the child grain      
                euler_id = find_closest_orientation_hcp(euler_id, euler_0);
                gPhi1_c(ind) = euler_id(1);
                gPhi_c(ind) = euler_id(2);
                gPhi2_c(ind) = euler_id(3);
                
                misorientation(jj,1) = calculate_misorientation_euler_d(euler_0, euler_id, 'hcp');
                col_ID = [col_ID; id];
                col_euler = [col_euler; euler_id];
            end
            tbl = table(col_ID, col_euler, misorientation); % for debugging  
            
            % (1) Find out which 'children' is actually the 'parent', by checking the misorientation 
            % Note: multiple 'grains' can have the parent orientation
            inds = find(misorientation < po_tolerance_angle);   
            if ~isempty(inds)
                % Use the average of all the 'child grains' with parent orientation.
                id_po = id_c(inds);
                inds_po = find(ismember(gID_c,id_po));
                euler_po = calculate_average_dominant_euler_hcp([gPhi1_c(inds_po), gPhi_c(inds_po), gPhi2_c(inds_po)]);
                    
                % remove the 'child grains' with the 'parent orientation', leave only the true 'twin children'
                id_c(inds) = [];
            else
                % the child maybe fully twinned. Just use euler_0 as euler_po  
                euler_po = euler_0;
                warning('The grain might have completely twinned');
                str = sprintf('iE = %d, ii = %d, ID_parent = %d \n', iE, ii, id_p);
                disp(str);
            end            
                
            % (2) Find out if the remaining child grain is a twin 
            for jj = 1:length(id_c)
                id = id_c(jj);  % twin grains id
                ind = (gID_c == id);
                euler_id = [gPhi1_c(ind), gPhi_c(ind), gPhi2_c(ind)];
                
                misorientation = [];
                %  if (demoTF==1)&&(jj==1)&&(ii==ii_demo)
                %       hcp_cell('euler',euler_id,'material','Mg','plotPlane',0,'plotBurgers',0,'plotTrace',0);
                %  end
                
                % Determine if this child grain is a twin variant:  
                % Compare [euler of this child] vs [euler of the iTwin system of the parent orientation]  
                for kk = 1:6
                    euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                    misorientation(kk) = calculate_misorientation_euler_d(euler_id, euler_kk, 'HCP');
                    % if (demoTF==1)&&(jj==1)&&(ii==ii_demo)
                    % 	hcp_cell('euler',euler_po,'material','Mg','ss',24+kk);
                    % 	hcp_cell('euler',euler_kk,'material','Mg','ss',24+kk,'plotPlane',0,'plotBurgers',0,'plotTrace',0);
                    % end
                end
                [min_val, iVariant_child] = min(abs(misorientation));    % find out
                
                % ==============> The child grain may be a twin area containing multiple variants. Assume the child orientation represents at least one true twin orientation  
                % Ff small enough, the child grain should be a twin. Do point-wise analysis   
                if min_val < twin_tolerance_angle
                    ID_variant_grain_wise(ID_c == id) = iVariant_child; % grain-wise variant map 
                    
                    ind_list = find(ID_c==id);
                    for kk = 1:length(ind_list)
                        ind = ind_list(kk);
                        euler_c = [phi1_c(ind), phi_c(ind), phi2_c(ind)]; 
                        
                        misorientation = [];
                        % compare [euler of this pixel]  vs  [euler of the iTwin system of the parent orientation]  
                        for kk = 1:6
                            euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                            misorientation(kk) = calculate_misorientation_euler_d(euler_c, euler_kk, 'HCP');
                        end
                        [miso, iVariant] = min(abs(misorientation));
                        Misorientation_point_wise(ind) = miso;
                        ID_variant_point_wise(ind) = iVariant;
                    end                    
                else
                    warning('Twin grain misorientation with parent orientation > 10 deg, rejected as a variant:');
                    str = sprintf('iE = %d, ii = %d, ID_parent = %d, jj = %d, ID_twin = %d \n', iE, ii, id_p, jj, id);
                    disp(str);
                end
            end
        end
    end
    variant_grain_wise{iE} = ID_variant_grain_wise;
    variant_point_wise{iE} = ID_variant_point_wise;
    
    % After processing this iE, Update the modified grain data and save
    save_dir_5 = [save_dir,'\step-5'];
    mkdir(save_dir_5);
    copyfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), ...
            fullfile(save_dir_5, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    gPhi1 = gPhi1_c;
    gPhi = gPhi_c;
    gPhi2 = gPhi2_c;
    save(fullfile(save_dir_5, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), 'gPhi1','gPhi','gPhi2', '-append');
    

    myplot(x,y, ID_variant_grain_wise, boundary_p); caxis([0 6]);
    title(['iE=',num2str(iE)]);
    print(fullfile(save_dir,['variant_grain_wise_iE=',num2str(iE),'.tif']),'-dtiff');
    close;
    
    myplot(x,y, ID_variant_point_wise, boundary_p); caxis([0 6]);
    title(['iE=',num2str(iE)]);
    print(fullfile(save_dir,['variant_pt_wise_iE=',num2str(iE),'.tif']),'-dtiff');
    close;      

end

save(fullfile(save_dir,'variant_maps.mat'),'variant_grain_wise','variant_point_wise');

% copy the updated file back to the main folder
for iE = 1:13
    copyfile(fullfile(save_dir_5, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), ...
        fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
end


%% part-6: summarize twin pct
load(fullfile(save_dir,'variant_maps.mat'),'variant_grain_wise','variant_point_wise');
twinPct = zeros(12,14);
for iE = 1:13
    load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), 'ID');

    variant_map = variant_point_wise{iE};
    if iE==1
        [nR,nC] = size(variant_map);
    end
    nr = 3;
    nc = 4;
    for ir=1:nr
       for ic = 1:nc
           sub_variant_map = variant_map([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic);
           sub_ID_map = ID([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic); % just for counting number of pixels
           
           ii = iE + 1;
           twinPct((ir-1)*nc+ic,ii) = sum(sub_variant_map(:)>0)/sum(sub_ID_map(:)>0);
       end
    end
end

tAvg = mean(twinPct);
tStd = std(twinPct);

save(fullfile(save_dir, 'twin_pct.mat'), 'twinPct', 'tAvg', 'tStd');

%% part-7, calculate EBSD estimated strain
load(fullfile(save_dir,'geotrans_and_id_link.mat'),'tforms');
for iE = 0:13
    iB = iE+1;
    if ~isempty(tforms{iB})
        [T,R,Z,S] = decompose_affine2d(tforms{iB});
        strain_ebsd(iB) = round(Z(1)-1, 4);
    else
        strain_ebsd(iB) = 0;
    end
end
str = [sprintf('strain_ebsd = ['), sprintf('%.4f, ',strain_ebsd(1:4)), newline, ...
    sprintf('%.4f, ',strain_ebsd(5:8)), newline, ...
    sprintf('%.4f, ',strain_ebsd(9:11)), newline, ...
    sprintf('%.4f, ',strain_ebsd(12:13)), sprintf('%.4f];',strain_ebsd(14))];
disp(str);

%% plot
load(fullfile(save_dir, 'twin_pct.mat'), 'twinPct', 'tAvg', 'tStd');
% strain for iE=0:13
% EBSD strain from fine alignment
strain_ebsd = [0, -0.0085, -0.0220, -0.0342, ...
    -0.0343, -0.0268, -0.0150, -0.0032, ...
    -0.0112, -0.0207, -0.0297, ...
    -0.0235, -0.0130, -0.0045];

strain_sg = [0, -0.0075, -0.015, -0.025, ...
    -0.023, -0.017, -0.0075, 0.002, ...
    -0.0075, -0.015, -0.025, ...
    -0.017, -0.0075, 0];

colors = parula(5);

% [1] using strain gage strain
figure; hold on;
inds = {1:3, 3:7, 7:11, 11:14};

errorbar(strain_sg(inds{1}), 100*tAvg(inds{1}), 100*tStd(inds{1}), '.-', 'color',colors(1,:), 'linewidth',1.5,'markersize',24);
errorbar(strain_sg(inds{2}), 100*tAvg(inds{2}), 100*tStd(inds{2}), '.-', 'color',colors(2,:), 'linewidth',1.5,'markersize',24);
errorbar(strain_sg(inds{3}), 100*tAvg(inds{3}), 100*tStd(inds{3}), '.-', 'color',colors(3,:), 'linewidth',1.5,'markersize',24);
errorbar(strain_sg(inds{4}), 100*tAvg(inds{4}), 100*tStd(inds{4}), '.-', 'color',colors(4,:), 'linewidth',1.5,'markersize',24);

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',[-0.04, 0.005],'ylim',[-2 50],'fontsize',18,'fontweight','normal');
xlabel('Strain from strain gage');
ylabel('Twin Area Percent (%)');
print(fullfile(save_dir,'twin_pct_vs_sg.tiff'),'-dtiff');

% [2] using (fine transformed) ebsd estimated strain
figure; hold on;
inds = {1:3, 3:7, 7:11, 11:14};

errorbar(strain_ebsd(inds{1}), 100*tAvg(inds{1}), 100*tStd(inds{1}), '.-', 'color',colors(1,:), 'linewidth',1.5,'markersize',24);
errorbar(strain_ebsd(inds{2}), 100*tAvg(inds{2}), 100*tStd(inds{2}), '.-', 'color',colors(2,:), 'linewidth',1.5,'markersize',24);
errorbar(strain_ebsd(inds{3}), 100*tAvg(inds{3}), 100*tStd(inds{3}), '.-', 'color',colors(3,:), 'linewidth',1.5,'markersize',24);
errorbar(strain_ebsd(inds{4}), 100*tAvg(inds{4}), 100*tStd(inds{4}), '.-', 'color',colors(4,:), 'linewidth',1.5,'markersize',24);

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',[-0.04, 0.005],'ylim',[-2 50],'fontsize',18,'fontweight','normal');
xlabel('Strain from ebsd estimate');
ylabel('Twin Area Percent (%)');
print(fullfile(save_dir,'twin_pct_vs_ebsd_strain.tiff'),'-dtiff');


tbl = array2table([(0:13)', strain_sg(:), strain_ebsd(:), 100*tAvg(:), 100*tStd(:)]);
tbl.Properties.VariableNames = {'iE','strain_sg','strain_ebsd','twinPct %','twinStd %'};
disp(tbl);
figure;
uitable('Data',tbl{:,:},'ColumnName',tbl.Properties.VariableNames,...
    'RowName',tbl.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
print(fullfile(save_dir,'twin pct table.tiff'),'-dtiff');

save(fullfile(save_dir, 'twin_pct.mat'), 'twinPct', 'tAvg', 'tStd', 'tbl');

%%
% % [3] using (rough transformed) ebsd estimate strain. Note, this sample is noisy, so the estimation might differ a lot ..   
% % EBSD strain from fine alignment
% strain_ebsd_0 = [0, -0.0079, -0.0213, -0.0335, ...
%     -0.0327, -0.0250, -0.0140, -0.0031, ...
%     -0.0111, -0.0202, -0.0286, ...
%     -0.0232, -0.0125, -0.0045];
% 
% figure; hold on;
% inds = {1:3, 3:7, 7:11, 11:14};
% 
% errorbar(strain_ebsd_0(inds{1}), 100*tAvg(inds{1}), 100*tStd(inds{1}), '.-', 'color',colors(1,:), 'linewidth',1.5,'markersize',24);
% errorbar(strain_ebsd_0(inds{2}), 100*tAvg(inds{2}), 100*tStd(inds{2}), '.-', 'color',colors(2,:), 'linewidth',1.5,'markersize',24);
% errorbar(strain_ebsd_0(inds{3}), 100*tAvg(inds{3}), 100*tStd(inds{3}), '.-', 'color',colors(3,:), 'linewidth',1.5,'markersize',24);
% errorbar(strain_ebsd_0(inds{4}), 100*tAvg(inds{4}), 100*tStd(inds{4}), '.-', 'color',colors(4,:), 'linewidth',1.5,'markersize',24);
% 
% set(gca,'xdir','normal','linewidth',1.5);
% set(gca,'xlim',[-0.04, 0.005],'ylim',[-2 50],'fontsize',18,'fontweight','normal');
% xlabel('Strain from ebsd estimate');
% ylabel('Twin Area Percent (%)');
% print(fullfile(save_dir,'twin_pct_vs_ebsd_0_strain.tiff'),'-dtiff');






