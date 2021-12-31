% analyze persistent twinning, using overlap


%% [1] Mg4Al_U2
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD';
save_dir = [working_dir, '\analysis'];
sample_name = 'Mg4Al_U2';

output_dir = ['E:\zhec umich Drive\0_temp_output\2021-05-04 Analyze persistent twin\', [sample_name, ' reMake']];
mkdir(output_dir);

% strain of this sample
strain_ebsd = [0.0000, -0.0085, -0.0218, -0.0336, ...
    -0.0337, -0.0265, -0.0149, -0.0032, ...
    -0.0111, -0.0205, -0.0293, ...
    -0.0232, -0.0130, -0.0045];

%% load ref data at iE = 0
iE = 0;
d = load(fullfile(save_dir, 'Mg4Al_U2_parent_grain_file_iE_0.mat'));
ID_0 = d.ID;
x = d.x;
y = d.y;
boundary_0 = find_one_boundary_from_ID_matrix(ID_0);

ID_overlap = ID_0;

%% load twin variant data
d = load(fullfile(save_dir, 'variant_maps.mat'));
variant_point_wise = d.variant_point_wise;
variant_grain_wise = d.variant_grain_wise;
for iE = 0:13
   iB = iE + 1;
   variant_pixel_iB_cell{iB} = variant_point_wise{iB};
   variant_grain_iB_cell{iB} = variant_grain_wise{iB};
end

% load the tforms
load(fullfile(save_dir, 'geotrans_and_id_link.mat'), 'tforms');
tforms{1} = fitgeotrans([0 0; 0 1; 1 0], [0 0; 0 1; 1 0], 'affine');  % assign at iE=0
%% load parent/child grain data at iEs, transform to iE=0
for iE = 0:13
    iB = iE + 1;
    
    % parent grain data
    d = load(fullfile(save_dir, ['Mg4Al_U2_parent_grain_file_iE_',num2str(iE),'.mat']));
    ID_p = d.ID;
    boundary_p = find_one_boundary_from_ID_matrix(ID_p);
    
    ID_p_iB_to_1_cell{iB} = interp_data(x,y,ID_p, x,y, tforms{iB}.invert, 'interp', 'nearest');
    boundary_p_iB_to_1_cell{iB} = interp_data(x,y,boundary_p, x,y, tforms{iB}.invert, 'interp', 'nearest');
    
    % child grain data
    d = load(fullfile(save_dir, ['Mg4Al_U2_grain_file_iE_',num2str(iE),'.mat']));
    ID_c = d.ID;
    boundary_c = find_one_boundary_from_ID_matrix(ID_c);
    
    ID_c_iB_to_1_cell{iB} = interp_data(x,y,ID_c, x,y, tforms{iB}.invert, 'interp', 'nearest');
    boundary_c_iB_to_1_cell{iB} = interp_data(x,y,boundary_c, x,y, tforms{iB}.invert, 'interp', 'nearest');
    
    % variant maps
    variant_pixel_iB_to_1_cell{iB} = interp_data(x,y,variant_pixel_iB_cell{iB}, x,y, tforms{iB}.invert, 'interp', 'nearest');
    variant_grain_iB_to_1_cell{iB} = interp_data(x,y,variant_grain_iB_cell{iB}, x,y, tforms{iB}.invert, 'interp', 'nearest');
    
    ID_overlap(ID_overlap~=ID_p_iB_to_1_cell{iB}) = 0;
end

gList = nan_unique(ID_overlap(:));  % gIDs in area of interest
gList(gList==0) = [];
gList(isnan(gList)) = [];

% check boundary map
for ii=[]
    figure; hold on;
    colors = parula(20);
    for iE = 0:13
        iB = iE + 1;
        ind = boundary_p_iB_to_1_cell{iB}==1;
        plot(x(ind), y(ind), '.', 'color',colors(iB,:));
    end
end

%% get variantID in geotransformed + valid area, for each iE
for iE = 0:13
    iB = iE + 1;
    map = variant_pixel_iB_to_1_cell{iB};
    map(ID_overlap==0) = 0;
    variant_pixel_cell{iB} = map;   % [map] variant number 1-6, transformed, valid area   
    current_twin_cell{iB} = double(variant_pixel_cell{iB}>0);  % [map] twin/not twin (type_double) 0 or 1, transformed, valid area  
    
    % cat variant_pixel in 3D to take mode. But first change variant=0 to
    % nan, to prevent 0 be counted by function mode.
    map(map==0) = nan;  
    if iE==0
        mapz = map;
    else
        mapz = cat(3,mapz,map);
    end
end
% the major variant identified at each pixel, if that pixel is ever twinned  
variant_pixel_mode = mode(mapz,3);
variant_pixel_mode(isnan(variant_pixel_mode)) = 0;  % change back to 0

%% (1) Pixel level summary, at each iE, there are:
% (p1) current twin (directly from twin map)
% (p2) past-or-present twin (obtained by accummulating twin map up to this iE) 
% based on these two basic maps, we can get
% (p3, yellow) fresh-twin = currently twinned - ever twinned before = p1 - p2(iE-1);  
% (p4, blue) recurring-twin = currently twinned & ever twinned bofore = p1 & p2(iE-1); 
% (p5, gray) de-twin = ever twinned before - currently twinned = p2(iE-1) - p1;
% p1 = p3 + p4
% p2 = p3 + p4 + p5
% 
% So, we have these pixel level maps:
% [1] p1 map, current twin map;
% [2] p2 map, ever twinned map; (p1 is a subset of p2)
% [3] p3 + p4 + p5 map, fre/recurring/de-twin map; 
p0_notwin = 0;
p1_current = 1;
p2_past_or_present = 2;
p3_fresh = 3;
p4_recur = 4;
p5_detwin = 5;

for iE = 0:13
    iB = iE + 1;
    if iE==0
        past_or_present_twin_cell{iB} = current_twin_cell{iB}>0;  % zeros(size(ID_0));    % pixels have twinned up to iE
        fresh_twin_cell{iB} = zeros(size(ID_0));    
        recurring_twin_cell{iB} = current_twin_cell{iB}>0;      % [Q1]===> existing polishing twins, as 'recurring twin' rather than 'fresh twin' pixel?
        de_twin_cell{iB} = zeros(size(ID_0));  
    else
        past_or_present_twin_cell{iB} = double(past_or_present_twin_cell{iB-1} | current_twin_cell{iB}>0);
        fresh_twin_cell{iB} = double(current_twin_cell{iB}>0 & past_or_present_twin_cell{iB-1}==0);
        recurring_twin_cell{iB} = double(current_twin_cell{iB}>0 & past_or_present_twin_cell{iB-1}>0);
        de_twin_cell{iB} = double(current_twin_cell{iB}==0 & past_or_present_twin_cell{iB-1}>0);
    end
end

% new/re/de-twin map: p3=new-twin, p4=re-twin, p5=de-twin
for iE = 0:13    
    iB = iE + 1;
    
    frd_twin_cell{iB} = fresh_twin_cell{iB}*p3_fresh + recurring_twin_cell{iB}*p4_recur + de_twin_cell{iB}*p5_detwin;
    
    % Plot to illustrate
    [f,a,c] = myplot(x,y,frd_twin_cell{iB}, boundary_p_iB_to_1_cell{iB});
        
    colors = parula(16);
    cmap = [1 1 1;      % backgrond, value < 2
        colors(14,:);   % yellowis, new twin, value = 3
        colors(2,:);    % blueish, retwin, value=4
        .5, .5, .5];    % gray, detwin, value=5
    caxis([1.5, 5.5]);
    colormap(cmap);
    
    set(gca,'xTickLabel',[],'yTickLabel',[]);
    set(c,'limits',[2.5,5.5], 'Ticks',[3,4,5], 'TickLabels',{'p3: Fresh Twin Pixel', 'p4: Recurring Twin Pixel', 'p5: Detwin Pixel'}); 
    title(['\fontsize{16}Load Step ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    set(gca,'xTickLabel',[],'yTickLabel',[], 'fontsize',16);
    set(gcf,'position', get(gcf,'position') .* [1 1 0 0] + [0 0 1000 600]);
    print(fullfile(output_dir,['frd twinned iE=',num2str(iE),'.tiff']),'-dtiff');
    
    close all;
end


%% (2) Grain level summary
% type-1 grain, current twin grain:
g1_new = 1; % (g1) new twin type-1 grain
g2_detwin_retwin = 2; % (g2) detwin then retwin type-1 grain
g3_evolving_1 = 3; % (g3) evolving type-1 grain
% type-2 grain, past-or-present twin grain (ever twinned up to iE)
g4_evolving_2 = 4; % (g4) evolving type-2 twin
g5_detwin_2 = 5; % (g5) completely detwin type-2 grain

[nR,nC] = size(ID_0);

for iE = 0:13
    iB = iE + 1;
    
    type_1_grain_label = zeros(nR,nC);    
    if iE==0
        type_2_grain_ID = zeros(nR,nC);
        type_2_grain_label = zeros(nR,nC);
    else
        type_2_grain_ID = type_2_grain_ID_cell{iB-1};
        type_2_grain_label = type_2_grain_label_cell{iB-1};
    end
    ID_assign = max(type_2_grain_ID(:)) + 1;  % initial ID to be assigned to ever twin grain 
    
    ID_c = ID_c_iB_to_1_cell{iB};  
    ID_c(ID_overlap==0) = 0;    % ONLY use valid/overlap area !
    
    % [step-1] loop check every 'current twin grain'
    twin_grain_list = nan_unique(ID_c(:));
    twin_grain_list(isnan(twin_grain_list)) = [];   % remove nan
    twin_grain_list(twin_grain_list==0) = [];       % remove 0 (if any)
    
    disp(['iE=',num2str(iE),', # of grains=',num2str(length(twin_grain_list))]);
    pause(1);
    for ii=1:length(twin_grain_list)
       ID_c_current = twin_grain_list(ii);
             
       % index of ID_c_current
       inds = ismember(ID_c, ID_c_current);
       
       % If this grain is a twin grain (contains twin pixel)
       grain_twinned_TF = any(current_twin_cell{iB}(inds));
       if grain_twinned_TF 
           if iE == 0
               type_1_grain_label(inds) = g3_evolving_1;    % [Q1]===> existing twin, evolving grain?
           else
               % Task: determine current twin map labels: g1 = new twin grain, g2 = detwin then retwinned, g3 = evolving current twin
               % to determine label, look at the pixel level label at previous load step
               labels = frd_twin_cell{iB-1}(inds);
               if sum(labels==p0_notwin)/numel(labels) >0.9    % intersect(labels, [p0,p3,p4,p5]) == p0
                   % the grain is completely new twin (g1/new twin), as it does not overlap with any previously twinned pixel
                   type_1_grain_label(inds) = g1_new;
               elseif intersect(labels, [p3_fresh, p4_recur, p5_detwin]) == p5_detwin
                   % if most are p5, we can allow a little bit p0 pixels
                   type_1_grain_label(inds) = g2_detwin_retwin;
               else
                   % evolving grain
                   type_1_grain_label(inds) = g3_evolving_1;
               end
           end
           
           % Task: generate past_or_present_twin_grain_ID: current twin grain will contribute to the 'ever twin grain map' 
           ids = unique(type_2_grain_ID(inds));
           ids(ids==0) = [];
           if isempty(ids)
               % if current twin grain does not overlap with any ever-twinned grain, ADD this grain to the ever twin grain map 
               type_2_grain_ID(inds) = ID_assign;
           else
               % if this grain overlap with ever twinned grain, MERGE them. (modify inds to include all grains)  
               inds = ismember(type_2_grain_ID, ids) | inds;
               type_2_grain_ID(inds) = ID_assign;
           end
           ID_assign = ID_assign + 1;
           
           % this 'merged ever twin grain' contains 'current twin', so at least p4/re-twin, maybe p3/new-twin pixels
           type_2_grain_label(inds) = g4_evolving_2;    % 'g4/evolving'
           
       end
    end
    type_1_grain_label_cell{iB} = type_1_grain_label;
    
    % [step-2] loop check every 'ever twin grain', find the completely detwinned grain ===> maybe need to move this to [step-1] 
    type_2_grain_list = nan_unique(type_2_grain_ID(:));
    type_2_grain_list(isnan(type_2_grain_list)) = []; % remove nan
    type_2_grain_list(type_2_grain_list==0) = [];     % remove 0 (if any)  
    
    for ii = 1:length(type_2_grain_list)
        ID_current = type_2_grain_list(ii);
        inds = ismember(type_2_grain_ID, ID_current);        
           
        % Task: make past_or_present_twin_grain_label   
        % Check the pixel level label at the current load step
        labels = unique(frd_twin_cell{iB}(inds));
        if intersect(labels, [p3_fresh,p4_recur,p5_detwin]) == p5_detwin
            % only contains conpletely detwinned pixels (p5), this ever twin grain should be labeled by (g4)   
            type_2_grain_label(inds) = g5_detwin_2;    % 'g5/completely detwin'
        else
            % this currently twinned contains p3/new-twin or p4/re-twin pixels.   
            type_2_grain_label(inds) = g4_evolving_2;    % 'g4/mixed'
        end
    end
    
    if iE>0
        % remake IDs
        type_2_grain_ID = hungarian_assign_ID_map(type_2_grain_ID, type_2_grain_ID_cell{iB-1});
    end
    type_2_grain_ID_cell{iB} = type_2_grain_ID;
    type_2_grain_label_cell{iB} = type_2_grain_label;
    
end

%% plot to check ==> OK.
close all;
for iE = 0:13
    iB = iE + 1;
    
    map = zeros(nR,nC);
    
    inds = type_2_grain_label_cell{iB}>0;
    map(inds) = type_2_grain_label_cell{iB}(inds);
    
    inds = type_1_grain_label_cell{iB}>0;
    map(inds) = type_1_grain_label_cell{iB}(inds);
    
    [f,a,c] = myplot(x,y, map, boundary_p_iB_to_1_cell{iB});
    
    colors = plasma(25);
    cmap = zeros(6,3);
    cmap(1,:) = [1,1,1];    % 0 = untwinned, background
    cmap(2,:) = colors(24,:);    % g1 = completely new twin, yellowish
    cmap(3,:) = colors(10,:);     % g2 = de-twin then re-twin, purple
    cmap(4,:) = colors(16,:);     % g3 = evoling current twin, pink
    cmap(5,:) = [0.7, 0.7, 0.7];     % g4 = evolving past-or-present twin, light gray
    cmap(6,:) = [0.2, 0.2, 0.2];    % g5 = completely de-twin, dark gray
    
    colormap(cmap);
    caxis([-0.5, 5.5]);
    set(c,'limits',[0.5, 5.5], 'Ticks', 1:5, ...
        'TickLabels', {'g1: New Twin Type-1 Grain','g2: Detwin Then Retwin Type-1 Grain','g3: Evolving Type-1 Grain', ...
            'g4: Evolving Type-2 Grain','g5: Completely Detwin Type-2 Grain'})
    title(['\fontsize{16}Load Step ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    set(gcf,'position',[80,172,1200,600]);
    set(gca,'xTickLabel',[],'yTickLabel',[],'fontsize',16);
    print(fullfile(output_dir,['grain label iE=',num2str(iE),'.tiff']),'-dtiff');
    close;
end

%% save data
save(fullfile(output_dir, 'Mg4Al_U2 twin data.mat'), 'gList', 'ID_overlap', ...
    'ID_p_iB_to_1_cell', 'ID_c_iB_to_1_cell', ...
    'boundary_p_iB_to_1_cell', 'boundary_c_iB_to_1_cell', ...
    'variant_grain_iB_to_1_cell', 'variant_pixel_iB_to_1_cell', ...
    'variant_pixel_cell', 'variant_pixel_mode', ...
    'frd_twin_cell', ...
    'current_twin_cell', 'type_1_grain_label_cell', ...
    'past_or_present_twin_cell', 'type_2_grain_ID_cell', 'type_2_grain_label_cell');


