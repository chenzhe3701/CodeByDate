%% setup
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-01-15 UM_134 Mg_C2 insitu EBSD';
save_dir = [working_dir, '\analysis'];
cd(working_dir);

%% reference, iE=0
iE = 0;
d = load(fullfile(save_dir, ['UM134_Mg_C2_grain_file_iE_',num2str(iE),'.mat']),'ID','phi1','phi','phi2','x','y',...
    'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');

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

%% Select deformed iE to select control grains [No need to run every time]
for iE = 1:13
close all;

% load the modified EBSD data
load(fullfile(save_dir, ['UM134_Mg_C2_modified_parent_grain_file_iE_',num2str(iE),'.mat']));

[boundary, boundaryID, neighborID, tripleTF, tripleID, indTriple, triIDs] = find_one_boundary_from_ID_matrix(ID);

% plot to help selection
myplot(x, y, ID_0, grow_boundary(boundary_0)); 
title('ID, iE=0','fontweight','normal');
set(gca,'fontsize',18);
label_map_with_ID(x,y,ID_0,gcf, [12, 13, 143,136],'r',24,10);  % illustrate selected control grains  
% print(fullfile(save_dir, ['b_ID_map_iE=0.tif']),'-r300','-dtiff');

myplot(x, y, ID, grow_boundary(boundary)); 
title(['ID, iE=',num2str(iE)],'fontweight','normal');
set(gca,'fontsize',18);
print(fullfile(save_dir, ['b_ID_map_iE=',num2str(iE),'.tif']),'-r300','-dtiff');
% label_map_with_ID(x,y,ID,gcf,[18,26,205, 71],'r',24,10);    % illustrate selected control grains

%% Double check if the grains for geotrans from iE=0 -> iE has changed.
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

% [[remove]] grains that are shown in BOTH the col_2 of id_link and col_2 of id_link_additional  
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
g_0 = id_link(:,1);
g_iE = id_link(:,2);  
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
% ===> just for illustration, I don't want to show the newly assigned ID
inds = ismember(ID_new, id_link_additional(:,1));
ID_new(inds) = nan;

myplot(x,y, ID_new, boundary_0_to_iE);
title(['iE=0 transformed with matched ID at iE=',num2str(iE)],'fontweight','normal');
set(gca,'fontsize',18);
print(fullfile(save_dir, ['b_tform_matched_ID_map_iE=',num2str(iE),'.tif']),'-r300','-dtiff');
% text(210,170,10000,'70<=>71','fontsize',18,'color','r')

end

%% Task 1: Affine transform ID map. Link grains.
% Generate geotrans/tform information at iE = 1 to 13, without showing results
% Link ids (save in 'tbl') at differnt iEs after loading the saved geotrans/tform. 
%   So, later, we can have grain files output again with different IDs.

for iE = 1:13
    disp(['iE=',num2str(iE)]);
    % load the modified EBSD data at iE
    load(fullfile(save_dir, ['UM134_Mg_C2_modified_parent_grain_file_iE_',num2str(iE),'.mat']));
    
    % match IDs
    g_0 = grain_pair{iE}(:,1);    % ref at iE=0
    g_iE = grain_pair{iE}(:,2);    % iE > 0, considered as deformed
    [~, loc_0] = ismember(g_0, gID_0);
    [~, loc_iE] = ismember(g_iE, gID);
    cpFrom = [gCenterX_0(loc_0),gCenterY_0(loc_0)];     % cpFrom is from ref image (iE=0)
    cpTo = [gCenterX(loc_iE),gCenterY(loc_iE)];         % cpTo is from deformed image (iE>0)
    
    % [1st] rough tform. The tform is to transform the [ref @ iE=0] to [deformed iE>0]
    tform = fitgeotrans(cpFrom, cpTo, 'affine');    % [x_0, y_0, 1] * tform.T = [x_iE, y_iE, 1]
    ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');

    [ID_new, id_link_additional, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);
    
    % [[remove]] '0's from linked ids
    ind_r = find(id_link(:,1)==0 | id_link(:,2)==0);
    id_link(ind_r,:) = [];
    inds = isnan(ID_0_to_iE)|(ID_new==0); % old ID is nan, or matched ID=0
    ID_new(inds) = nan;

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
    tforms{iE} = tform;     % save the tform
    
    ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');
   
    [ID_new, id_link_additional, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);
    % remove 0s
    ind_r = find(id_link(:,1)==0 | id_link(:,2)==0);
    id_link(ind_r,:) = [];
    inds = isnan(ID_0_to_iE)|(ID_new==0); % old ID is nan, or matched ID=0
    ID_new(inds) = nan;

    tbl_iE = array2table(id_link);
    tbl_iE.Properties.VariableNames = {'iE_0',['iE_',num2str(iE)]};
    
    if iE==1
        tbl = tbl_iE;
    else
        tbl = innerjoin(tbl,tbl_iE);
    end
end

save(fullfile(save_dir,'geotrans_and_id_link.mat'),'tforms','tbl');


%% Task 2: Modify Euler angles, find closest orientation, and [change ID] to the matched IE_0
% For all the matched grains at different iEs, change grain euler angles to
% the crystallographically-equivalent one that has the closest orientation
% to that of the ref grain when calculated without considering symmetry.

load(fullfile(save_dir,'geotrans_and_id_link.mat'),'tforms','tbl');

maxID = 0;
for iE = 0:13
    if iE==0
        f = matfile(fullfile(save_dir,['UM134_Mg_C2_grain_file_iE_',num2str(iE)]));
    else
        f = matfile(fullfile(save_dir,['UM134_Mg_C2_modified_parent_grain_file_iE_',num2str(iE)]));
    end
    gID = f.gID;
    maxID = max(maxID, max(gID));
end
id_to_add = 10^(ceil(log10(maxID)));    % round to 1000 etc for adjustment. 

for iE = 0:13
    clear euler_aligned_to_sample;
    if iE==0
        load(fullfile(save_dir,['UM134_Mg_C2_grain_file_iE_',num2str(iE)]));
    else        
        load(fullfile(save_dir,['UM134_Mg_C2_modified_parent_grain_file_iE_',num2str(iE)]));
    end
    if ~exist('euler_aligned_to_sample','var')
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
        
        eulerd = find_closest_orientation_hcp(eulerd, eulerd_ref);      % ==> Find Closest Orientation
        
        gPhi1(gID==id_source) = eulerd(1);
        gPhi(gID==id_source) = eulerd(2);
        gPhi2(gID==id_source) = eulerd(3);
        
        gID(gID==id_source) = id_target;
        gNeighbors(gNeighbors==id_source) = id_target;
    end
    
    boundary = find_boundary_from_ID_matrix(ID);
    myplot(x,y,ID,boundary);
    clim = caxis;
    caxis([0,maxID]);
    title(['iE=',num2str(iE),' matched ID #']);
    print(fullfile(save_dir, ['matched_ID_map_iE=',num2str(iE),'.tif']),'-r300','-dtiff');
    close;
    
    euler_aligned_to_parent = 1;
    save(fullfile(save_dir, ['UM134_Mg_C2_parent_grain_file_iE_',num2str(iE),'.mat']), 'euler_aligned_to_sample','euler_aligned_to_parent', ...
        'ID','phi1','phi','phi2','x','y',...
        'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');

end

