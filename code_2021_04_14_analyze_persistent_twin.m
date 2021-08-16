% analyze persistent twinning, using overlap


%% [1] Mg4Al_U2
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD';
save_dir = [working_dir, '\analysis'];
sample_name = 'Mg4Al_U2';

output_dir = ['E:\zhec umich Drive\0_temp_output\2021-04-20 Analyze persistent twin\',sample_name];
mkdir(output_dir);

%% [2] Mg4Al_C1
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD';
save_dir = [working_dir, '\analysis'];
sample_name = 'Mg4Al_C3';

output_dir = ['E:\zhec umich Drive\0_temp_output\2021-04-20 Analyze persistent twin\',sample_name];
mkdir(output_dir);

%% [3] UM134_C1, OK but not good
% clear; clc; close all;
% addChenFunction;
% working_dir = 'E:\zhec umich Drive\2020-12-05 UM134 Mg_C1 insitu EBSD';
% save_dir = [working_dir, '\analysis'];
% sample_name = 'UM134_Mg_C1';
% 
% output_dir = ['E:\zhec umich Drive\0_temp_output\2021-04-20 Analyze persistent twin\',sample_name];
% mkdir(output_dir);

%% [3] UM134_C2
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-01-15 UM134 Mg_C2 insitu EBSD';
save_dir = [working_dir, '\analysis'];
sample_name = 'UM134_Mg_C2';

output_dir = ['E:\zhec umich Drive\0_temp_output\2021-04-20 Analyze persistent twin\',sample_name];
mkdir(output_dir);

%% [4] UM134_C3
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-01-29 UM134 Mg_C3 insitu EBSD'
save_dir = [working_dir, '\analysis'];
sample_name = 'UM134_Mg_C3';

output_dir = ['E:\zhec umich Drive\0_temp_output\2021-04-20 Analyze persistent twin\',sample_name];
mkdir(output_dir);



%% load variant maps, tforms
% variant maps as cell data
d = matfile(fullfile(save_dir, 'variant_maps.mat'));
variant_map_cell = d.variant_point_wise;

% tform data
d = matfile(fullfile(save_dir,'geotrans_and_id_link'));
tforms = d.tforms;

%% load reference data at iE = 0, reference ID map and grain boundary
iE = 0;
d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
ID_0 = d.ID;
gPhi1_0 = d.gPhi1;
gPhi_0 = d.gPhi;
gPhi2_0 = d.gPhi2;
x = d.x;
y = d.y;

boundary_0 = find_one_boundary_from_ID_matrix(ID_0);

%%
for iE = 1:13
    iB = iE + 1;
    
    % data at iE
    d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    ID = d.ID;
    boundary = find_one_boundary_from_ID_matrix(ID);
    
    tform = tforms{iB}; % tform is from i0 to iE
    variant_iE = variant_map_cell{iE};
    
    % transform maps at iE to i0
    variant_iE_to_0 = interp_data(x,y,variant_iE, x,y, tform.invert, 'interp', 'nearest');
    boundary_iE_to_0 = interp_data(x,y,boundary, x,y, tform.invert, 'interp', 'nearest');
    ID_iE_to_0 = interp_data(x,y,ID, x,y, tform.invert, 'interp', 'nearest');
    
    % record the maps in cell
    variant_iE_to_0_cell{iE} = variant_iE_to_0;
    boundary_iE_to_0_cell{iE} = boundary_iE_to_0;
    ID_iE_to_0_cell{iE} = ID_iE_to_0;
    
%     boundary_cell{iE} = boundary;
%     ID_cell{iE} = ID;
end


%%
% for each grain, calculate SchmidFactor of each variant
% For each pixel, (1) when does it first twin, (2) when does it first
% de-twin, (3) when does it first re-twin, (4) when does it de-twin again

% unique ID list, i.e., IDs that exist in all iEs
gList = unique(ID_0(:));
for iE = 1:13
    ID = ID_iE_to_0_cell{iE};
    gList = intersect(gList, unique(ID(:)));
end

% ID_overlap, where IDs at all iEs are the same without overlap error due to geometric transform  
ID_overlap = ID_0;
ID_overlap(~ismember(ID_overlap, gList)) = 0;
for iE = 1:13
    ID_overlap(ID_overlap ~= ID_iE_to_0_cell{iE}) = 0;
end

% valid variant map, valid in areas where ID overlap in all iEs
for iE = 1:13
    temp_map = variant_iE_to_0_cell{iE};
    temp_map(ID_overlap==0) = 0;
    variant_cell{iE} = temp_map;
    twin_cell{iE} = double(temp_map>0);
end

% initialize
max_twin_cell{1} = double(twin_cell{1}>0);    % max-twin area up to iE
new_twin_cell{1} = double(twin_cell{1}>0);    % new-twin at iE. max_twin(iE-1) + new_twin(iE) = max_twin(iE)  
re_twin_cell{1} = zeros(size(ID));  % re_twin(iE) + new_twinned(iE) = twin(iE)  
de_twin_cell{1} = zeros(size(ID));  % max_twin(iE) - de_twin(iE) = twin(iE)

for iE = 2:13
    max_twin_cell{iE} = double(max_twin_cell{iE-1} | twin_cell{iE}>0);
    new_twin_cell{iE} = double(twin_cell{iE}>0 & max_twin_cell{iE-1}==0);
    re_twin_cell{iE} = double(twin_cell{iE}>0 & max_twin_cell{iE-1}>0);
    de_twin_cell{iE} = double(max_twin_cell{iE}>0 & twin_cell{iE}==0);
end


%% just for illustration (1) max-twin, (2) new-twin, (3) re-twin, (4) de-twin area
for iE = 1:13
    close all;
    myplot(max_twin_cell{iE}, boundary_iE_to_0_cell{iE}); title(['max-twin iE=',num2str(iE)], 'fontweight','normal');
    print(fullfile(output_dir,['max-twin iE=',num2str(iE),'.tiff']),'-dtiff');
    
    myplot(new_twin_cell{iE}, boundary_iE_to_0_cell{iE}); title(['new-twin iE=',num2str(iE)], 'fontweight','normal');
    print(fullfile(output_dir,['new-twin iE=',num2str(iE),'.tiff']),'-dtiff');
    
    myplot(re_twin_cell{iE}, boundary_iE_to_0_cell{iE}); title(['re-twin iE=',num2str(iE)], 'fontweight','normal');
    print(fullfile(output_dir,['re-twin iE=',num2str(iE),'.tiff']),'-dtiff');
    
    myplot(de_twin_cell{iE}, boundary_iE_to_0_cell{iE}); title(['de-twin iE=',num2str(iE)], 'fontweight','normal');
    print(fullfile(output_dir,['de-twin iE=',num2str(iE),'.tiff']),'-dtiff');
end

close all;


%% This might be the only useful plot
for iE = 1:13
    % myplot(twin_cell{iE} - re_twin_cell{iE} - new_twin_cell{iE});   % this confirms: new_twin + re_twin == twin
    close all;
    temp_map = zeros(size(ID));
    inds = re_twin_cell{iE}==1;
    temp_map(inds) = temp_map(inds) + 1;
    
    inds = new_twin_cell{iE}==1;
    temp_map(inds) = temp_map(inds) + 2;
    
    inds = de_twin_cell{iE}==1;
    temp_map(inds) = temp_map(inds) + 3;
    
    [f,a,c] = myplot(temp_map, boundary_iE_to_0_cell{iE});
    
    colors = parula(16);
    cmap = [1,1,1; 
        colors(2,:);
        colors(14,:);
        .5, .5, .5]; % back ground row 1, de-twin color last row
    
    colormap(cmap);
    caxis([-0.5, 3.5]);
    set(c,'limits',[0.5,3.5], 'Ticks',[1,2,3], 'TickLabels',{'re-twin', 'new-twin', 'de-twin'});
    % set(c,'limits',[0.5,2.5], 'Ticks',[1,2], 'TickLabels',{sprintf('%s\\newline%s', 're-twin', '(or retained)'),'new-twin'});
    title(['iE=',num2str(iE)],'fontweight','norma');
    print(fullfile(output_dir,['re and new twinned iE=',num2str(iE),'.tiff']),'-dtiff');
end

%%
% for ii = 1:length(gList)
%     ID_current = gList(ii);
%     % [indR_min, indR_max, indC_min, indC_max] = find_inds_local(ID, ID_current);
%     
%     ID = ID_iE_to_0_cell{iE};
%     temp_map = variant_iE_to_0_cell{iE};
%     inds = ismember(ID, ID_current);  % inds of current grain at iE
%     variants = temp_map(inds);
% end









