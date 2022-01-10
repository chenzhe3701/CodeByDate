%% Generate pixel-level and grain-level twin category summary for in-situ EBSD data
clear; clc; close all;
addChenFunction;

check_boundary_overlay_quality = 1;

% [data_dir, sample_name, sample_ID, plot_symbol, group_number, ymax]
cells = {'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1, 300;
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1, 300;
    'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2, 350;
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2, 350;
    'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3, 200;
    'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3, 200;
    % 'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'x', 4, 450;
    'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 'o', 4, 450;
    'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 's', 4, 450;
    % 'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'x', 5, 200;
    'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 'o', 5, 200;
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 's', 5, 200};

variant_map_dir = 'E:\zhec umich Drive\All twin variant maps cleaned';


% Where to save the output figures
output_dir = 'E:\zhec umich Drive\0_temp_output\all twin evolution analysis';
mkdir(output_dir);

%% Make twin type maps, pixel level summary and grain level summary
for icell = 1:size(cells,1)
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    sample_ID = cells{icell,3};
    
    % load ref data at iE = 0
    d = load(fullfile(sample_dir, [sample_name, '_parent_grain_file_iE_0.mat']));
    ID_0 = d.ID;
    x = d.x;
    y = d.y;
    boundary_0 = find_one_boundary_from_ID_matrix(ID_0);
    
    ID_overlap = ID_0;
    
    % load twin variant data
    d = load(fullfile(variant_map_dir, [sample_name,'_variant_maps.mat']));
    variant_point_wise = d.variant_point_wise;
    variant_grain_wise = d.variant_grain_wise;
    for iE = 0:13
        iB = iE + 1;
        vp_iB_cell{iB} = variant_point_wise{iB};
        vg_iB_cell{iB} = variant_grain_wise{iB};
    end
    
    % load the tforms
    load(fullfile(sample_dir, 'geotrans_and_id_link.mat'), 'tforms');
    tforms{1} = fitgeotrans([0 0; 0 1; 1 0], [0 0; 0 1; 1 0], 'affine');  % assign at iE=0
    
    %% load parent/child grain data at iEs, transform to iE=0
    for iE = 0:13
        iB = iE + 1;
        
        % parent grain data
        d = load(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        boundary_p = find_one_boundary_from_ID_matrix(ID_p);
        
        ID_p_iB_to_1_cell{iB} = interp_data(x,y,ID_p, x,y, tforms{iB}.invert, 'interp', 'nearest');
        boundary_p_iB_to_1_cell{iB} = interp_data(x,y,boundary_p, x,y, tforms{iB}.invert, 'interp', 'nearest');
        
        % child grain data
        d = load(fullfile(sample_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));
        ID_c = d.ID;
        boundary_c = find_one_boundary_from_ID_matrix(ID_c);
        
        ID_c_iB_to_1_cell{iB} = interp_data(x,y,ID_c, x,y, tforms{iB}.invert, 'interp', 'nearest');
        boundary_c_iB_to_1_cell{iB} = interp_data(x,y,boundary_c, x,y, tforms{iB}.invert, 'interp', 'nearest');
        
        % variant maps
        vp_iB_to_1_cell{iB} = interp_data(x,y,vp_iB_cell{iB}, x,y, tforms{iB}.invert, 'interp', 'nearest');
        vg_iB_to_1_cell{iB} = interp_data(x,y,vg_iB_cell{iB}, x,y, tforms{iB}.invert, 'interp', 'nearest');
        
        ID_overlap(ID_overlap~=ID_p_iB_to_1_cell{iB}) = 0;
    end
    
    gList = nan_unique(ID_overlap(:));  % gIDs in area of interest
    gList(gList==0) = [];
    gList(isnan(gList)) = [];
    
    % check boundary map
    if check_boundary_overlay_quality
        figure; hold on;
        colors = parula(20);
        for iE = 0:13
            iB = iE + 1;
            ind = boundary_p_iB_to_1_cell{iB}==1;
            plot(x(ind), y(ind), '.', 'color',colors(iB,:));
        end
        axis equal;
        print(fullfile(output_dir,[sample_name,' gb_oly_quality.tiff']),'-dtiff');
        close;
    end
    
    %% get variantID in geotransformed + valid area, for each iE
    for iE = 0:13
        iB = iE + 1;

        % I think we need to further clean up the map. Because non-valid areas are removed, there might be some small isolated 'twin grains' left.  
        % purpose: in overlap area, connected twin area of a child grain should be large enough  
        map_t = ID_c_iB_to_1_cell{iB};   % map_t = all child grains
        twinTF = vp_iB_to_1_cell{iB}>0; % twin area
        map_t(~twinTF) = 0;  % now, map_t = twinned child grains  
        map_t(ID_overlap==0) = 0;   % now, map_t = twinned child grains in overlap area 
        map_t = one_pass_clean(map_t,16);   % now, map_t = twinned child grains in overlap area, and are big enough  
        big_twin_grain_in_overlap_area_cell{iB} = map_t;
        
        vp_map = vp_iB_to_1_cell{iB};
        vp_map(big_twin_grain_in_overlap_area_cell{iB}==0) = 0;
        
        variant_pixel_cell{iB} = vp_map;   % [map] variant number 1-6, transformed, valid area
        current_twin_cell{iB} = double(variant_pixel_cell{iB}>0);  % [map] twin/not twin (type_double) 0 or 1, transformed, valid area
        
        % cat variant_pixel in 3D to take mode. But first change variant=0 to
        % nan, to prevent 0 be counted by function mode.
        vp_map(vp_map==0) = nan;
        if iE==0
            mapz = vp_map;
        else
            mapz = cat(3,mapz,vp_map);
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
        % title(['\fontsize{16}Load Step ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','normal');
        title(['\fontsize{16}Load Step ',num2str(iE)],'fontweight','normal');
        set(gca,'xTickLabel',[],'yTickLabel',[], 'fontsize',16);
        set(gcf,'position', get(gcf,'position') .* [1 1 0 0] + [0 0 1000 600]);
        print(fullfile(output_dir,[sample_name ' frd twinned iE=',num2str(iE),'.tiff']),'-dtiff');
        
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
        
        type_1_grain_ID = zeros(nR,nC);
        type_1_grain_label = zeros(nR,nC);
        if iE==0
            type_2_grain_ID = zeros(nR,nC);
            type_2_grain_label = zeros(nR,nC);
        else
            type_2_grain_ID = type_2_grain_ID_cell{iB-1};
            type_2_grain_label = type_2_grain_label_cell{iB-1};
        end
        ID2_assign = max(type_2_grain_ID(:)) + 1;  % initial ID to be assigned to type-2 grain
        
        ID_c = ID_c_iB_to_1_cell{iB};
        ID_c(big_twin_grain_in_overlap_area_cell{iB}==0) = 0;    % ONLY use valid/overlap area !
        
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
                type_1_grain_ID(inds) = ID_c_current;
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
                
                % Task: generate type-2 grain ID: current twin grain will contribute to the type-2 grain map
                ids = unique(type_2_grain_ID(inds));
                ids(ids==0) = [];
                if isempty(ids)
                    % if current twin grain does not overlap with any type-2 grain, ADD this grain to the ever twin grain map
                    type_2_grain_ID(inds) = ID2_assign;
                else
                    % if this grain overlap with type-2 grain, MERGE them. (modify inds to include all grains)
                    inds = ismember(type_2_grain_ID, ids) | inds;
                    type_2_grain_ID(inds) = ID2_assign;
                end
                ID2_assign = ID2_assign + 1;
                
                % this 'merged type-2 grain' contains 'current twin', so at least p4/re-twin, maybe p3/new-twin pixels
                type_2_grain_label(inds) = g4_evolving_2;    % 'g4/evolving'
                
            end
        end
        type_1_grain_ID_cell{iB} = type_1_grain_ID;
        type_1_grain_label_cell{iB} = type_1_grain_label;
        
        % [step-2] loop check every type-2 grain, find the completely detwinned grain ===> maybe need to move this to [step-1]
        type_2_grain_list = nan_unique(type_2_grain_ID(:));
        type_2_grain_list(isnan(type_2_grain_list)) = []; % remove nan
        type_2_grain_list(type_2_grain_list==0) = [];     % remove 0 (if any)
        
        for ii = 1:length(type_2_grain_list)
            ID_current = type_2_grain_list(ii);
            inds = ismember(type_2_grain_ID, ID_current);
            
            % Task: make type-2 grain label
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
        % title(['\fontsize{16}Load Step ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','normal');
        title(['\fontsize{16}Load Step ',num2str(iE)],'fontweight','normal')
        set(gcf,'position',[80,172,1200,600]);
        set(gca,'xTickLabel',[],'yTickLabel',[],'fontsize',16);
        print(fullfile(output_dir,[sample_name, ' grain label iE=',num2str(iE),'.tiff']),'-dtiff');
        close;
    end
    
    %% save data
    save(fullfile(output_dir, [sample_name, ' twin evolution label.mat']), 'gList', 'ID_overlap', 'big_twin_grain_in_overlap_area_cell', ...
        'ID_p_iB_to_1_cell', 'ID_c_iB_to_1_cell', ...
        'boundary_p_iB_to_1_cell', 'boundary_c_iB_to_1_cell', ...
        'vg_iB_to_1_cell', 'vp_iB_to_1_cell', ...
        'variant_pixel_cell', 'variant_pixel_mode', ...
        'current_twin_cell', 'past_or_present_twin_cell', 'frd_twin_cell', ...
        'type_1_grain_ID_cell', 'type_1_grain_label_cell', ...
        'type_2_grain_ID_cell', 'type_2_grain_label_cell');
    
end