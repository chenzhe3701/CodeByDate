% analyze persistent twinning, using overlap
working_dir = 'E:\zhec umich Drive\0_temp_output\2021-04-20 Analyze persistent twin';

%% [1] Mg4Al_C1
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD';
save_dir = [working_dir, '\analysis'];
cd(save_dir);

sample_name = 'Mg4Al_C1';

% variant maps as cell data
d = matfile(fullfile(save_dir, 'variant_maps.mat'));
variant_map_cell = d.variant_point_wise;

% tform data
d = matfile(fullfile(save_dir,'geotrans_and_id_link'));
tforms = d.tforms;

% iE = 0, reference ID map and grain boundary
iE = 0;
d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
ID_0 = d.ID;
x = d.x;
y = d.y;

boundary_0 = find_one_boundary_from_ID_matrix(ID_0);

% make and save inverse transformed maps at all iEs
for iE = 1:7
    iB = iE + 1;
    
    % data at iE
    d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    ID = d.ID;
    boundary = find_one_boundary_from_ID_matrix(ID);
    
    tform = tforms{iB}; % tform is from i0 to iE
    vMap_iE = variant_map_cell{iE};
    
    % transform maps at iE to i0
    vMap_iE_to_0 = interp_data(x,y,vMap_iE, x,y, tform.invert, 'interp', 'nearest');
    boundary_iE_to_0 = interp_data(x,y,boundary, x,y, tform.invert, 'interp', 'nearest');
    ID_iE_to_0 = interp_data(x,y,ID, x,y, tform.invert, 'interp', 'nearest');
    
    % record the maps in cell
    vMap_iE_to_0_cell{iE} = vMap_iE_to_0;
    boundary_iE_to_0_cell{iE} = boundary_iE_to_0;
    ID_iE_to_0_cell{iE} = ID_iE_to_0;
    
    boundary_cell{iE} = boundary;
    ID_cell{iE} = ID;
end

% [[[check]]] the quality of transform. Overlay all iEs
close all;
figure; hold on;
z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');

colors = linspecer(7);
for iE = 1:7
  z = boundary_iE_to_0_cell{iE};
  z(z==0) = nan;
  plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color',colors(iE,:));
end
axis equal;
set(gca,'xlim', [min(x(:)),max(x(:))], 'ylim',[min(y(:)),max(y(:))]);
title([strrep(sample_name,'_',' '),' gb all iEs to iE=0, overlay'], 'fontweight','normal');

% [[[check]]] deformed iE=3 overlay with iE=0
figure; hold on;
z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');
z = boundary_cell{3};
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'r');
axis equal;
set(gca,'xlim', [min(x(:)),max(x(:))], 'ylim',[min(y(:)),max(y(:))]);
title([strrep(sample_name,'_',' '), ' gb iE=0(black) and iE=3(red)'], 'fontweight','normal');

%% [2] Mg4Al_C3
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD';
save_dir = [working_dir, '\analysis'];
cd(save_dir);

sample_name = 'Mg4Al_C3';

% variant maps as cell data
d = matfile(fullfile(save_dir, 'variant_maps.mat'));
variant_map_cell = d.variant_point_wise;

% tform data
d = matfile(fullfile(save_dir,'geotrans_and_id_link'));
tforms = d.tforms;

% iE = 0, reference ID map and grain boundary
iE = 0;
d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
ID_0 = d.ID;
x = d.x;
y = d.y;

boundary_0 = find_one_boundary_from_ID_matrix(ID_0);

% make and save inverse transformed maps at all iEs
for iE = 1:13
    iB = iE + 1;
    
    % data at iE
    d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    ID = d.ID;
    boundary = find_one_boundary_from_ID_matrix(ID);
    
    tform = tforms{iB}; % tform is from i0 to iE
    vMap_iE = variant_map_cell{iE};
    
    % transform maps at iE to i0
    vMap_iE_to_0 = interp_data(x,y,vMap_iE, x,y, tform.invert, 'interp', 'nearest');
    boundary_iE_to_0 = interp_data(x,y,boundary, x,y, tform.invert, 'interp', 'nearest');
    ID_iE_to_0 = interp_data(x,y,ID, x,y, tform.invert, 'interp', 'nearest');
    
    % record the maps in cell
    vMap_iE_to_0_cell{iE} = vMap_iE_to_0;
    boundary_iE_to_0_cell{iE} = boundary_iE_to_0;
    ID_iE_to_0_cell{iE} = ID_iE_to_0;
    
    boundary_cell{iE} = boundary;
    ID_cell{iE} = ID;
end

% [[[check]]] the quality of transform. Overlay all iEs
close all;
figure; hold on;
z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');

colors = linspecer(13);
for iE = 1:13
  z = boundary_iE_to_0_cell{iE};
  z(z==0) = nan;
  plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color',colors(iE,:));
end
axis equal;
set(gca,'xlim', [min(x(:)),max(x(:))], 'ylim',[min(y(:)),max(y(:))]);
title([strrep(sample_name,'_',' '),' gb all iEs to iE=0, overlay'], 'fontweight','normal');

% [[[check]]] deformed iE=3 overlay with iE=0
figure; hold on;
z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');
z = boundary_cell{3};
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'r');
axis equal;
set(gca,'xlim', [min(x(:)),max(x(:))], 'ylim',[min(y(:)),max(y(:))]);
title([strrep(sample_name,'_',' '), ' gb iE=0(black) and iE=3(red)'], 'fontweight','normal');


%% [3] Mg4Al_U2
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD';
save_dir = [working_dir, '\analysis'];
cd(save_dir);

sample_name = 'Mg4Al_U2';

% variant maps as cell data
d = matfile(fullfile(save_dir, 'variant_maps.mat'));
variant_map_cell = d.variant_point_wise;

% tform data
d = matfile(fullfile(save_dir,'geotrans_and_id_link'));
tforms = d.tforms;

% iE = 0, reference ID map and grain boundary
iE = 0;
d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
ID_0 = d.ID;
x = d.x;
y = d.y;

boundary_0 = find_one_boundary_from_ID_matrix(ID_0);

% make and save inverse transformed maps at all iEs
for iE = 1:13
    iB = iE + 1;
    
    % data at iE
    d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    ID = d.ID;
    boundary = find_one_boundary_from_ID_matrix(ID);
    
    tform = tforms{iB}; % tform is from i0 to iE
    vMap_iE = variant_map_cell{iE};
    
    % transform maps at iE to i0
    vMap_iE_to_0 = interp_data(x,y,vMap_iE, x,y, tform.invert, 'interp', 'nearest');
    boundary_iE_to_0 = interp_data(x,y,boundary, x,y, tform.invert, 'interp', 'nearest');
    ID_iE_to_0 = interp_data(x,y,ID, x,y, tform.invert, 'interp', 'nearest');
    
    % record the maps in cell
    vMap_iE_to_0_cell{iE} = vMap_iE_to_0;
    boundary_iE_to_0_cell{iE} = boundary_iE_to_0;
    ID_iE_to_0_cell{iE} = ID_iE_to_0;
    
    boundary_cell{iE} = boundary;
    ID_cell{iE} = ID;
end

% [[[check]]] the quality of transform. Overlay all iEs
close all;
figure; hold on;
z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');

colors = linspecer(13);
for iE = 1:13
  z = boundary_iE_to_0_cell{iE};
  z(z==0) = nan;
  plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color',colors(iE,:));
end
axis equal;
set(gca,'xlim', [min(x(:)),max(x(:))], 'ylim',[min(y(:)),max(y(:))]);
title([strrep(sample_name,'_',' '),' gb all iEs to iE=0, overlay'], 'fontweight','normal');

% [[[check]]] deformed iE=3 overlay with iE=0
figure; hold on;
z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');
z = boundary_cell{3};
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'r');
axis equal;
set(gca,'xlim', [min(x(:)),max(x(:))], 'ylim',[min(y(:)),max(y(:))]);
title([strrep(sample_name,'_',' '), ' gb iE=0(black) and iE=3(red)'], 'fontweight','normal');



%% [4] UM134_C1
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2020-12-05 UM134 Mg_C1 insitu EBSD';
save_dir = [working_dir, '\analysis'];
cd(save_dir);

sample_name = 'UM134_Mg_C1';

% variant maps as cell data
d = matfile(fullfile(save_dir, 'variant_maps.mat'));
variant_map_cell = d.variant_point_wise;

% tform data
d = matfile(fullfile(save_dir,'geotrans_and_id_link'));
tforms = d.tforms;

% iE = 0, reference ID map and grain boundary
iE = 0;
d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
ID_0 = d.ID;
x = d.x;
y = d.y;

boundary_0 = find_one_boundary_from_ID_matrix(ID_0);

% make and save inverse transformed maps at all iEs
for iE = 1:13
    iB = iE + 1;
    
    % data at iE
    d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    ID = d.ID;
    boundary = find_one_boundary_from_ID_matrix(ID);
    
    tform = tforms{iB}; % tform is from i0 to iE
    vMap_iE = variant_map_cell{iE};
    
    % transform maps at iE to i0
    vMap_iE_to_0 = interp_data(x,y,vMap_iE, x,y, tform.invert, 'interp', 'nearest');
    boundary_iE_to_0 = interp_data(x,y,boundary, x,y, tform.invert, 'interp', 'nearest');
    ID_iE_to_0 = interp_data(x,y,ID, x,y, tform.invert, 'interp', 'nearest');
    
    % record the maps in cell
    vMap_iE_to_0_cell{iE} = vMap_iE_to_0;
    boundary_iE_to_0_cell{iE} = boundary_iE_to_0;
    ID_iE_to_0_cell{iE} = ID_iE_to_0;
    
    boundary_cell{iE} = boundary;
    ID_cell{iE} = ID;
end

% [[[check]]] the quality of transform. Overlay all iEs
close all;
figure; hold on;
z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');

colors = linspecer(13);
for iE = 1:13
  z = boundary_iE_to_0_cell{iE};
  z(z==0) = nan;
  plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color',colors(iE,:));
end
axis equal;
set(gca,'xlim', [min(x(:)),max(x(:))], 'ylim',[min(y(:)),max(y(:))]);
title([strrep(sample_name,'_',' '),' gb all iEs to iE=0, overlay'], 'fontweight','normal');

% [[[check]]] deformed iE=3 overlay with iE=0
figure; hold on;
z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');
z = boundary_cell{3};
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'r');
axis equal;
set(gca,'xlim', [min(x(:)),max(x(:))], 'ylim',[min(y(:)),max(y(:))]);
title([strrep(sample_name,'_',' '), ' gb iE=0(black) and iE=3(red)'], 'fontweight','normal');

%% [5] UM134_C2
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-01-15 UM134 Mg_C2 insitu EBSD';
save_dir = [working_dir, '\analysis'];
cd(save_dir);

sample_name = 'UM134_Mg_C2';

% variant maps as cell data
d = matfile(fullfile(save_dir, 'variant_maps.mat'));
variant_map_cell = d.variant_point_wise;

% tform data
d = matfile(fullfile(save_dir,'geotrans_and_id_link'));
tforms = d.tforms;

% iE = 0, reference ID map and grain boundary
iE = 0;
d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
ID_0 = d.ID;
x = d.x;
y = d.y;

boundary_0 = find_one_boundary_from_ID_matrix(ID_0);

% make and save inverse transformed maps at all iEs
for iE = 1:13
    iB = iE + 1;
    
    % data at iE
    d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    ID = d.ID;
    boundary = find_one_boundary_from_ID_matrix(ID);
    
    tform = tforms{iB}; % tform is from i0 to iE
    vMap_iE = variant_map_cell{iE};
    
    % transform maps at iE to i0
    vMap_iE_to_0 = interp_data(x,y,vMap_iE, x,y, tform.invert, 'interp', 'nearest');
    boundary_iE_to_0 = interp_data(x,y,boundary, x,y, tform.invert, 'interp', 'nearest');
    ID_iE_to_0 = interp_data(x,y,ID, x,y, tform.invert, 'interp', 'nearest');
    
    % record the maps in cell
    vMap_iE_to_0_cell{iE} = vMap_iE_to_0;
    boundary_iE_to_0_cell{iE} = boundary_iE_to_0;
    ID_iE_to_0_cell{iE} = ID_iE_to_0;
    
    boundary_cell{iE} = boundary;
    ID_cell{iE} = ID;
end

% [[[check]]] the quality of transform. Overlay all iEs
close all;
figure; hold on;
z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');

colors = linspecer(13);
for iE = 1:13
  z = boundary_iE_to_0_cell{iE};
  z(z==0) = nan;
  plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color',colors(iE,:));
end
axis equal;
set(gca,'xlim', [min(x(:)),max(x(:))], 'ylim',[min(y(:)),max(y(:))]);
title([strrep(sample_name,'_',' '),' gb all iEs to iE=0, overlay'], 'fontweight','normal');

% [[[check]]] deformed iE=3 overlay with iE=0
figure; hold on;
z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');
z = boundary_cell{3};
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'r');
axis equal;
set(gca,'xlim', [min(x(:)),max(x(:))], 'ylim',[min(y(:)),max(y(:))]);
title([strrep(sample_name,'_',' '), ' gb iE=0(black) and iE=3(red)'], 'fontweight','normal');


%% [6] UM134_C3
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-01-29 UM134 Mg_C3 insitu EBSD';
save_dir = [working_dir, '\analysis'];
cd(save_dir);

sample_name = 'UM134_Mg_C3';

% variant maps as cell data
d = matfile(fullfile(save_dir, 'variant_maps.mat'));
variant_map_cell = d.variant_point_wise;

% tform data
d = matfile(fullfile(save_dir,'geotrans_and_id_link'));
tforms = d.tforms;

% iE = 0, reference ID map and grain boundary
iE = 0;
d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
ID_0 = d.ID;
x = d.x;
y = d.y;

boundary_0 = find_one_boundary_from_ID_matrix(ID_0);

% make and save inverse transformed maps at all iEs
for iE = 1:13
    iB = iE + 1;
    
    % data at iE
    d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    ID = d.ID;
    boundary = find_one_boundary_from_ID_matrix(ID);
    
    tform = tforms{iB}; % tform is from i0 to iE
    vMap_iE = variant_map_cell{iE};
    
    % transform maps at iE to i0
    vMap_iE_to_0 = interp_data(x,y,vMap_iE, x,y, tform.invert, 'interp', 'nearest');
    boundary_iE_to_0 = interp_data(x,y,boundary, x,y, tform.invert, 'interp', 'nearest');
    ID_iE_to_0 = interp_data(x,y,ID, x,y, tform.invert, 'interp', 'nearest');
    
    % record the maps in cell
    vMap_iE_to_0_cell{iE} = vMap_iE_to_0;
    boundary_iE_to_0_cell{iE} = boundary_iE_to_0;
    ID_iE_to_0_cell{iE} = ID_iE_to_0;
    
    boundary_cell{iE} = boundary;
    ID_cell{iE} = ID;
end

% [[[check]]] the quality of transform. Overlay all iEs
close all;
figure; hold on;
z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');

colors = linspecer(13);
for iE = 1:13
  z = boundary_iE_to_0_cell{iE};
  z(z==0) = nan;
  plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color',colors(iE,:));
end
axis equal;
set(gca,'xlim', [min(x(:)),max(x(:))], 'ylim',[min(y(:)),max(y(:))]);
title([strrep(sample_name,'_',' '),' gb all iEs to iE=0, overlay'], 'fontweight','normal');

% [[[check]]] deformed iE=3 overlay with iE=0
figure; hold on;
z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');
z = boundary_cell{3};
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'r');
axis equal;
set(gca,'xlim', [min(x(:)),max(x(:))], 'ylim',[min(y(:)),max(y(:))]);
title([strrep(sample_name,'_',' '), ' gb iE=0(black) and iE=3(red)'], 'fontweight','normal');


