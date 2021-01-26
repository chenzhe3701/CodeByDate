%% EBSD data for in-situ test, Mg4Al_C3, tested on 2020-12-23 (14 load steps, 2 cycles)
% Tasks:
% (1) Add boundary to selected grains and generate new EBSD/grain data
% (2) Align Euler angles to sample reference frame.

%% setup
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD';
save_dir = [working_dir, '\analysis'];
cd(working_dir);

%% reference, iE=0
iE = 0;
[EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['Mg4Al_C3 grain_file_type_1 iE=',num2str(iE),'.txt']));
[EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['Mg4Al_C3 grain_file_type_2 iE=',num2str(iE),'.txt']));

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

[EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['Mg4Al_C3 parent_grain_file_type_1 iE=',num2str(iE),'.txt']));
[EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['Mg4Al_C3 parent_grain_file_type_2 iE=',num2str(iE),'.txt']));

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
label_map_with_ID(x,y,ID_0,gcf, [12, 13, 143,136],'r',24,10);  % illustrate selected control grains  
print(fullfile(save_dir, ['a_ID_map_iE=0.tif']),'-r300','-dtiff');

myplot(x, y, ID, grow_boundary(boundary)); 
title(['ID, iE=',num2str(iE)],'fontweight','normal');
set(gca,'fontsize',18);
print(fullfile(save_dir, ['a_ID_map_iE=',num2str(iE),'.tif']),'-r300','-dtiff');
% label_map_with_ID(x,y,ID,gcf,[18,26,205, 71],'r',24,10);    % illustrate selected control grains

%% Geotrans the iE=0 map to overlay on iE>0 map, to overlay the same grains [no need to run every time]
% (step-1) Select no less than 3 grain pairs for control points, and record
grain_pair{1} = [12, 13;
    13, 14;
    143, 139;
    136, 134];
grain_pair{2} = [12, 12;
    13, 13;
    143, 143;
    136, 138];
grain_pair{3} = [12, 12;
    13, 15;
    143, 146;
    136, 141];
grain_pair{4} = [12, 12;
    13, 15;
    143, 145;
    136, 139];
grain_pair{5} = [12, 12;
    13, 14;
    143, 145;
    136, 140];
grain_pair{6} = [12, 12;
    13, 13;
    143, 137;
    136, 132];
grain_pair{7} = [12, 12;
    13, 13;
    143, 139;
    136, 134];
grain_pair{8} = [12, 12;
    13, 13;
    143, 138;
    136, 133];
grain_pair{9} = [12, 12;
    13, 13;
    143, 135;
    136, 130];
grain_pair{10} = [12, 12;
    13, 15;
    143, 142;
    136, 136];
grain_pair{11} = [12, 12;
    13, 13;
    143, 136;
    136, 131];
grain_pair{12} = [12, 13;
    13, 14;
    143, 163;
    136, 157];
grain_pair{13} = [12, 12;
    13, 13;
    143, 136;
    136, 131];



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



%% There might not be equal number of grains on the reference vs deformed maps. 
% (1) If there are more grains on the reference map than the deformed map,
% the extra grain on the reference map will not be able to find a match. If
% necessary, we want to select the grain on the deformed map, alter the
% tolerance angle to divide that grain into two (or more) grains. (2) If
% there are more grains on the deformed map, it's more difficult to check
% on the result 'deformed map with matched ID', but in this situation, the
% reference map might need to be processed, i.e., select a grain and reduce
% the tolerance angle to divide it into two grains.
% 
% Ref can merge twin, but might not be useful at this point
%% Record (1) the ID_list (on map_iE) to correct for each iE, and (2) the tolerance for each grain
ID_list{1} = [57,103,117,129];
tolerance_cell{1} = [5,5,5,5];
ID_list{2} = [56,71,86,104,122,129,133];
tolerance_cell{2} = [5,5,5,5,5,5,5];
ID_list{3} = [60,68,77,91,109,126,132,132,133,136];
tolerance_cell{3} = [5,5,5,5,4,5,5,4,5,5];
ID_list{4} = [61,66,75,89,106,123,130,134];
tolerance_cell{4} = [5,5,5,5,5,5,5,5];
ID_list{5} = [61,75,91,107,124,130,131,135];
tolerance_cell{5} = [5,5,5,5,5,3,5,5];
ID_list{6} = [32,54,66,99,115,127];
tolerance_cell{6} = [5,5,5,5,5,5];
ID_list{7} = [31,54,117,129];
tolerance_cell{7} = [5,5,5,5];
ID_list{8} = [32,53,116,128];
tolerance_cell{8} = [5,5,5,5];
ID_list{9} = [30,51,66,98,114,121,125];
tolerance_cell{9} = [5,5,5,5,5,5,5,5,5];
ID_list{10} = [57,63,72,87,93,105,120,127,131];
tolerance_cell{10} = [5,5,4,5,5,5,5,5,5];
ID_list{11} = [57,65,82,87,100,114,126];
tolerance_cell{11} = [5,5,5,5,5,5,5];
ID_list{12} = [32,56,114,134,150];
tolerance_cell{12} = [5,5,5,5,5]; 
ID_list{13} = [29,114,126];
tolerance_cell{13} = [5,5,5];
% load the saved masks
try
    load(fullfile(save_dir,['Mg4Al_C3_mask_cell.mat']),'mask_cell','ID_updated_cell');
catch
   mask_cell = cell(1,13); 
end
boundary_new = boundary;    
iN = 1;
ID_updated = ID;
next_gID = max(gID) + 1;

% In the next part, we need to check if the added boundary segment just
% create one addtional grain.

%% So, here, the function is to :
% (1) select the grain, (2) calculate local misorientation map, (3) find a
% line with the largest misorientation to divide the grain into two.

N = length(ID_list{iE});
while iN <= N
    try
        close(hf);
    end
    str = sprintf('iE=%d, iN=%d', iE, iN);
    disp(str);
    % example for debug, iE=1, grain 57
    ID_current = ID_list{iE}(iN);
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
    if (length(mask_cell{iE})>=iN) && ~isempty(mask_cell{iE}{iN})
        mask = mask_cell{iE}{iN};
    else
        myplot(boundary_local);
        myplotm(misorientation_max);
        caxis([tolerance_cell{iE}(iN), 100]);
        h = drawpolygon;
        customWait(h);
        mask = h.createMask();
        mask_cell{iE}{iN} = mask;
    end
    
    tolerance = tolerance_cell{iE}(iN); % tolerance for this grain, default to 5
    boundary_local_new = (misorientation_max > tolerance)&(mask==1);
    boundary_local_new = double(ID_local~=ID_current | boundary_local_new);
    
    ID_local_new = find_ID_map_from_boundary_map(boundary_local_new);
    if length(unique(ID_local_new(:))) ~= 2
        % redraw
        mask_cell{iE}{iN} = [];     % delete mask
        warning(' More than one additional number of grains is created, you might want to do it again!');
    else
        % If split into g1 and g2, let g1 = 0, g2 = next_gID - ID_current, other grains = 0
        ID_local_new(ID_local_new==1) = 0;
        ID_local_new(ID_local_new==2) = next_gID - ID_current;
        ID_local_new(ID_local~=ID_current) = 0;
        ID_local_new = ID_local_new + ID_local;     % after adding, g2 will be next_gID
        ID_updated(indR_min:indR_max, indC_min:indC_max) = ID_local_new;    % update ID_new
        
        next_gID = next_gID + 1;    % update next_gID and iN
        iN = iN + 1;
        
        % update boundary, for illustration purpose.
        boundary_local_new = (misorientation_max > tolerance)&(mask==1);
        boundary_local_new = double(boundary_local | boundary_local_new);
        boundary_new(indR_min:indR_max, indC_min:indC_max) = boundary_local_new;
        close;close;
        hf = myplot(boundary_local_new);
    end
    
    ID_updated_cell{iE} = ID_updated;
    save(fullfile(save_dir, ['Mg4Al_C3_mask_cell.mat']),'mask_cell','ID_updated_cell');
end


close;close;
myplot(boundary_new);

%% Task-1: Use ID_temp to add grain boundaries, summarize and save the EBSD data, for iE
% ID_updated: modified from the added grain boundaries
% ID: the target ID map at this iE

umPerDp = 1;    % micron per data point in EBSD data
data_in.umPerDp = umPerDp;
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

save(fullfile(save_dir, ['Mg4Al_C3_modified_parent_grain_file_iE_',num2str(iE),'.mat']),'ID','phi1','phi','phi2','x','y',...
    'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');

close all;

disp(['finished iE=',num2str(iE)]);

end
%% Task-2: Align Euler angles to sample reference frame, for both [parent grain file] and [regular grain file]. 
close all;

% regular grain file
for iE = 0:13    
    
    [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['Mg4Al_C3 grain_file_type_1 iE=',num2str(iE),'.txt']));
    [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['Mg4Al_C3 grain_file_type_2 iE=',num2str(iE),'.txt']));
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
    save(fullfile(save_dir, ['Mg4Al_C3_grain_file_iE_',num2str(iE),'.mat']),'euler_aligned_to_sample','ID','phi1','phi','phi2','x','y',...
        'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');
    
end

% modified for parent grain file
for iE = 1:13    
    clear euler_aligned_to_sample;
    load(fullfile(save_dir, ['Mg4Al_C3_modified_parent_grain_file_iE_',num2str(iE),'.mat']));
    
    if exist('euler_aligned_to_sample')&&(euler_aligned_to_sample==1)
        error('no need to align euler to sample again. check.');
    else
        [phi1, phi, phi2] = align_euler_to_sample(phi1, phi, phi2, 'non-mtex', 90,180,0);
        [gPhi1, gPhi, gPhi2] = align_euler_to_sample(gPhi1, gPhi, gPhi2, 'non-mtex', 90,180,0);
        euler_aligned_to_sample = 1;
    end
    
    
    save(fullfile(save_dir, ['Mg4Al_C3_modified_parent_grain_file_iE_',num2str(iE),'.mat']),'euler_aligned_to_sample','ID','phi1','phi','phi2','x','y',...
        'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');
    
end






