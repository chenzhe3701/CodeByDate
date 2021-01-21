%% EBSD data for in-situ test, UM134 Pure Mg_C1, tested on 2020-12-05 (14 load steps, 2 cycles)

%% setup
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2020-12-05 UM_134 Mg_C1 insitu EBSD';
save_dir = [working_dir, '\analysis'];
cd(working_dir);

%% reference, iE=0
iE = 0;
[EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_1.txt']));
[EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_2.txt']));

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
boundary_0 = find_boundary_from_ID_matrix(ID_0);

%% Select deformed iE to select control grains [No need to run every time]
iE = 13;
close all;
try
    [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_1_parent.txt']));
    [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_2_parent.txt']));
catch
    [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_1.txt']));
    [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_2.txt']));
end
gID = EBSD_data_2(:,column_index_2(1));
gCenterX = EBSD_data_2(:,column_index_2(5));
gCenterY = EBSD_data_2(:,column_index_2(6));
gEdge = EBSD_data_2(:,column_index_2(10));
ID = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
boundary = find_boundary_from_ID_matrix(ID);

% plot to help selection
myplot(x, y, ID_0, boundary_0); 
title('ID, iE=0','fontweight','normal');
set(gca,'fontsize',18);
label_map_with_ID(x,y,ID_0,gcf, [42,34, 267,271],'r',24,10);  % illustrate selected control grains  

myplot(x, y, ID, boundary); 
title(['ID, iE=',num2str(iE)],'fontweight','normal');
set(gca,'fontsize',18);
% label_map_with_ID(x,y,ID,gcf,[18,26,205, 71],'r',24,10);    % illustrate selected control grains

%% Geotrans the iE=0 map to overlay on iE>0 map, to overlay the same grains [no need to run every time]
% (step-1) Select no less than 3 grain pairs, and record
grain_pair{1} = [42, 33;
    34, 30;
    267, 252;
    271, 268];
grain_pair{2} = [42, 32;
    34, 33;
    267, 253;
    271, 273];
grain_pair{3} = [42, 37;
    34, 38;
    267, 261;
    271, 280];
grain_pair{4} = [42, 36;
    34, 33;
    267, 253;
    271, 264];
grain_pair{5} = [42, 40;
    34, 36;
    267, 255;
    271, 266];
grain_pair{6} = [42, 45;
    34, 36;
    267, 257;
    271, 263];
grain_pair{7} = [42, 44;
    34, 39;
    267, 256;
    271, 266];
grain_pair{8} = [42, 37;
    34, 34;
    267, 248;
    271, 263];
grain_pair{9} = [42, 29;
    34, 32;
    267, 241;
    271, 259];
grain_pair{10} = [42, 26;
    34, 31;
    267, 239;
    271, 262];
grain_pair{11} = [42, 33;
    34, 34;
    267, 243;
    271, 260];
grain_pair{12} = [42, 37;
    34, 34;
    267, 243;
    271, 260];
grain_pair{13} = [42, 42;
    34, 36;
    267, 251;
    271, 261];

% (step-2) Use the selected control grains to roughly align. The result is already decent 
g_0 = grain_pair{iE}(:,1);    % ref at iE=0
g_iE = grain_pair{iE}(:,2);    % iE > 0, considered as deformed

[~, g_ind_0] = ismember(g_0, gID_0);
[~, g_ind_iE] = ismember(g_iE, gID);
cpFrom = [gCenterX_0(g_0),gCenterY_0(g_0)];   % cpFrom is from ref image (iE=0)
cpTo = [gCenterX(g_iE),gCenterY(g_iE)]; % cpTo is from deformed image (iE>0)

% The tform is to transform the [ref @ iE=0] to [deformed iE>0]
tform = fitgeotrans(cpFrom, cpTo, 'affine');    % [x_0, y_0, 1] * tform.T = [x_iE, y_iE, 1]
tform = maketform('affine', tform.T);   % make it into ftorm struct, so 'interp_data' can use
ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');
boundary_0_to_iE = find_boundary_from_ID_matrix(ID_0_to_iE);

myplot(x,y, ID_0_to_iE, boundary_0_to_iE); 
title(['(rough) affine transformed iE=0 -> iE=',num2str(iE)]);

[ID_new, id_link_additional, id_link_zero_overlap, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);
% remove '0's from linked ids
ind_r = find(id_link(:,1)==0 | id_link(:,2)==0);
id_link(ind_r,:) = [];


% (step-3) Based on id_link, use all non-edge grains to fit again, and show result! 
g_0 = id_link(:,1);
g_iE = id_link(:,2);    
% look at 'gEdge' property of grains in id_link, and remove edge grains in id_link 
ind = (gEdge_0(g_0)==1) | (gEdge(g_iE)==1);
g_0(ind) = [];
g_iE(ind) = [];

[~, g_ind_0] = ismember(g_0, gID_0);
[~, g_ind_iE] = ismember(g_iE, gID);

cpFrom = [gCenterX_0(g_0),gCenterY_0(g_0)]; % cpFrom is from deformed image (iE>0)
cpTo = [gCenterX(g_iE),gCenterY(g_iE)];   % cpTo is from ref image (iE=0)

% The tform is to transform the 'deformed' back into the 'ref'
tform = fitgeotrans(cpFrom, cpTo, 'affine');    % [x_ref, y_ref, 1] * tform.T = [x_def, y_def, 1]
tform = maketform('affine', tform.T);   % make it into ftorm struct, so 'interp_data' can use

ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');
boundary_0_to_iE = find_boundary_from_ID_matrix(ID_0_to_iE);

myplot(x,y, ID_0_to_iE, boundary_0_to_iE); 
title(['(final) affine transformed iE=0 -> iE=',num2str(iE)]);

[ID_new, id_link_additional, id_link_zero_overlap, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);
% remove 0s
ind_r = find(id_link(:,1)==0 | id_link(:,2)==0);
id_link(ind_r,:) = [];

myplot(x,y, ID_new, boundary_0_to_iE);
title(['iE=0 transformed with matched ID at iE=',num2str(iE)],'fontweight','normal');
set(gca,'fontsize',18);
% text(210,170,10000,'70<=>71','fontsize',18,'color','r')

% myplot(x,y, ID_new - ID, boundary);
% title(['ID\_affined\_reassigned - ID\_iE=',num2str(iE)]);

%% Generate geotrans information:  without showing results, just performe the link of all iEs from 1 to 7
for iE = 1:13
    % load iE
    try
        [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_1_parent.txt']));
        [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_2_parent.txt']));
    catch
        [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_1.txt']));
        [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_2.txt']));
    end
    gID = EBSD_data_2(:,column_index_2(1));
    gCenterX = EBSD_data_2(:,column_index_2(5));
    gCenterY = EBSD_data_2(:,column_index_2(6));
    gEdge = EBSD_data_2(:,column_index_2(10));
    ID = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
    boundary = find_boundary_from_ID_matrix(ID);
    
    % match IDs
    g_0 = grain_pair{iE}(:,1);    % ref at iE=0
    g_iE = grain_pair{iE}(:,2);    % iE > 0, considered as deformed
    
    [~, g_ind_0] = ismember(g_0, gID_0);
    [~, g_ind_iE] = ismember(g_iE, gID);
    cpFrom = [gCenterX_0(g_0),gCenterY_0(g_0)];   % cpFrom is from ref image (iE=0)
    cpTo = [gCenterX(g_iE),gCenterY(g_iE)]; % cpTo is from deformed image (iE>0)
    
    % The tform is to transform the 'deformed' back into the 'ref'
    tform = fitgeotrans(cpFrom, cpTo, 'affine');    % [x_ref, y_ref, 1] * tform.T = [x_def, y_def, 1]
    tform = maketform('affine', tform.T);   % make it into ftorm struct, so 'interp_data' can use
    ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');
    boundary_0_to_iE = find_boundary_from_ID_matrix(ID_0_to_iE);
    
    [ID_new, id_link_additional, id_link_zero_overlap, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);
    % remove 0s
    ind_r = find(id_link(:,1)==0 | id_link(:,2)==0);
    id_link(ind_r,:) = [];

    % (2) Based on id_link, use all non-edge grains to fit again
    g_0 = id_link(:,1);
    g_iE = id_link(:,2);
    % look at 'gEdge' property of grains in id_link, and remove edge grains in id_link
    ind = (gEdge_0(g_0)==1) | (gEdge(g_iE)==1);
    g_0(ind) = [];
    g_iE(ind) = [];
    
    [~, g_ind_0] = ismember(g_0, gID_0);
    [~, g_ind_iE] = ismember(g_iE, gID);
    
    cpFrom = [gCenterX_0(g_0),gCenterY_0(g_0)]; % cpFrom is from deformed image (iE>0)
    cpTo = [gCenterX(g_iE),gCenterY(g_iE)];   % cpTo is from ref image (iE=0)
    
    % The tform is to transform the 'deformed' back into the 'ref'
    tform = fitgeotrans(cpFrom, cpTo, 'affine');    % [x_ref, y_ref, 1] * tform.T = [x_def, y_def, 1]
    tform = maketform('affine', tform.T);   % make it into ftorm struct, so 'interp_data' can use
    
    % save the tform
    tforms{iE} = tform;
    
    ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');
    boundary_0_to_iE = find_boundary_from_ID_matrix(ID_0_to_iE);
    
    [ID_new, id_link_additional, id_link_zero_overlap, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);
    % remove 0s
    ind_r = find(id_link(:,1)==0 | id_link(:,2)==0);
    id_link(ind_r,:) = [];

    tbl_iE = array2table(id_link);
    tbl_iE.Properties.VariableNames = {'iE_0',['iE_',num2str(iE)]};
    
    if iE==1
        tbl = tbl_iE;
    else
        tbl = innerjoin(tbl,tbl_iE);
    end
end

save(fullfile(save_dir,'geotrans_and_id_link.mat'),'tforms');

%% Link ids at differnt iEs. Load save geotrans/tform, so we can have grain files output again with different IDs.
% load tforms
load(fullfile(save_dir,'geotrans_and_id_link.mat'),'tforms')

for iE = 1:13
    % load iE
    try
        [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_1_parent.txt']));
        [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_2_parent.txt']));
    catch
        [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_1.txt']));
        [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_2.txt']));
    end
    gID = EBSD_data_2(:,column_index_2(1));
    gCenterX = EBSD_data_2(:,column_index_2(5));
    gCenterY = EBSD_data_2(:,column_index_2(6));
    gEdge = EBSD_data_2(:,column_index_2(10));
    ID = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
    boundary = find_boundary_from_ID_matrix(ID);
    
    tform = tforms{iE};
    
    ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');
    boundary_0_to_iE = find_boundary_from_ID_matrix(ID_0_to_iE);
    
    [ID_new, id_link_additional, id_link_zero_overlap, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);
    % remove 0s
    ind_r = find(id_link(:,1)==0 | id_link(:,2)==0);
    id_link(ind_r,:) = [];

    tbl_iE = array2table(id_link);
    tbl_iE.Properties.VariableNames = {'iE_0',['iE_',num2str(iE)]};
    
    if iE==1
        tbl = tbl_iE;
    else
        tbl = innerjoin(tbl,tbl_iE);
    end
end

save(fullfile(save_dir,'geotrans_and_id_link.mat'),'tforms','tbl','-append');

%% Next step is to change grain euler angles

load(fullfile(save_dir,'geotrans_and_id_link.mat'),'tforms','tbl');

maxID = 0;
for iE = 0:13
    [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_2.txt']));
    gID = EBSD_data_2(:,column_index_2(1));
    maxID = max(maxID, max(gID));
end
ida = 10^(ceil(log10(maxID)));    % round to 1000 etc for adjustment

for iE = 0:13
    try
        [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_1_parent.txt']));
        [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_2_parent.txt']));
    catch
        [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_1.txt']));
        [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_2.txt']));
    end
    
    ID = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
    
    gID = EBSD_data_2(:,column_index_2(1));
    gPhi1 = EBSD_data_2(:,column_index_2(2));
    gPhi = EBSD_data_2(:,column_index_2(3));
    gPhi2 = EBSD_data_2(:,column_index_2(4));

    gNNeighbors = EBSD_data_2(:,column_index_2(7));
    gNeighbors = EBSD_data_2(:,(column_index_2(7)+1):(size(EBSD_data_2,2)));
    
    % adjust, add a big enough value to prevent old unmodified grains having the same id as new modified grains  
    ID(ID>0) = ID(ID>0) + ida;
    gID(gID>0) = gID(gID>0) + ida;
    gNeighbors(gNeighbors>0) = gNeighbors(gNeighbors>0) + ida;    % add
    
    id_link = tbl.Variables;
    g_list_target = id_link(:,1);
    g_list_source = id_link(:,iE+1) + ida;  % remember to use 'iE+1'
    for ii = 1:length(g_list_target)
        id_source = g_list_source(ii);
        id_target = g_list_target(ii);
        % modify
        ID(ID==id_source) = id_target;
        
        eulerd_ref = [gPhi1_0(gID_0==id_target), gPhi_0(gID_0==id_target), gPhi2_0(gID_0==id_target)];
        eulerd = [gPhi1(gID==id_source), gPhi(gID==id_source), gPhi2(gID==id_source)];
        eulerd = find_closest_orientation_hcp(eulerd, eulerd_ref);
        gPhi1(gID==id_source) = eulerd(1);
        gPhi(gID==id_source) = eulerd(2);
        gPhi2(gID==id_source) = eulerd(3);
        
        gID(gID==id_source) = id_target;
        gNeighbors(gNeighbors==id_source) = id_target;
    end
    
    boundary = find_boundary_from_ID_matrix(ID);
    myplot(x,y,ID,boundary);
    clim = caxis;
    caxis(clim-[0,ida]);
    title(['iE=',num2str(iE),' matched ID #']);
    
    save(fullfile(save_dir, ['data_with_ID_overwrite_iE_',num2str(iE),'.mat']), 'ID','gID','gPhi1','gPhi','gPhi2','gNeighbors');
end