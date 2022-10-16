% analyze detwin

clear; clc; close all;
addChenFunction;

% [data_dir, sample_name, sample_ID, plot_symbol, group_number, ymax, sample_material, ymax_for_detwin_variant_count, sample_label]
cells = {'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1, 50, 'Mg4Al FG', 160, 'Mg4Al fine grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1, 50, 'Mg4Al FG', 160, 'Mg4Al fine grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2, 350, 'Mg4Al MG', 160, 'Mg4Al medium grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2, 350, 'Mg4Al MG', 160, 'Mg4Al medium grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3, 100, 'Mg4Al CG', 80, 'Mg4Al coarse grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3, 100, 'Mg4Al CG', 80, 'Mg4Al coarse grain';
    % 'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'x', 4, 30, 'Mg FG', 80, 'Mg fine grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 'o', 4, 30, 'Mg FG', 80, 'Mg fine grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 's', 4, 30, 'Mg FG', 80, 'Mg fine grain';
    % 'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'x', 5, 200, 'Mg CG', 60, 'Mg coarse grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 'o', 5, 200, 'Mg CG', 60, 'Mg coarse grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 's', 5, 200, 'Mg CG', 60, 'Mg coarse grain'};

% variant_map_dir = 'H:\Other computers\My Laptop w541 2022\zhec umich Drive\All twin variant maps cleaned';

% location of the twin evolution (twin detwin retwin) data
evolution_dir = 'C:\Users\chenz\Work\Data\0_temp_output\all twin evolution analysis';

output_dir = 'C:\Users\chenz\Work\Data\0_temp_output\all twin evolution analysis';
mkdir(output_dir);

%% [analyze] detwin: variant pxls, max variant pxls, SF. Variant significantly detwinned or not.
% [Task-1] algorithm: to count pixel of each variant
% For each parent grain,
%  for each variant, max area = # pixels in type-2 grains
%   current area = # pixels in type-1 grains
%
% [Task-2] This is actually, pct of activated variants in the completely-detwinned type-2 grains.  
%   --> this might not be a very good analysis.
% algorithm: to determine 'significantly detwinned variant'
% For a completely detwinned type-2 grain, if variant pixel area >10% of type-2 grain area, then this is a 'significantly detwinned' variant.     
% Use reconstructed variant maps
% at iE, examine type-2 grain:
%  if completely-detwinned type-2 grain, determine variant, label as significantly detwinned variant, and was active at max compression
%  if evoling type-2 grain, determine variant, label as active variant at max compression
%   for each grain, have two 6x1 vectors, iv_active_at_max_compression, iv_detwinned_at_this_iE



mkdir(fullfile(output_dir, 'detwin'));

for icell = [2,5,7,10] % 1:size(cells,1)
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    sample_ID = cells{icell,3};
    sample_material = cells{icell,7};
    sample_label = cells{icell, 9};

    % load ref data at iE = 0, for SF calculation
    d = matfile(fullfile(sample_dir, [sample_name, '_parent_grain_file_iE_0.mat']));
    ID = d.ID;
    gPhi1 = d.gPhi1;
    gPhi = d.gPhi;
    gPhi2 = d.gPhi2;
    gID = d.gID;
    
    % twin evolution data
    d = matfile(fullfile(output_dir, [sample_name, ' twin evolution label.mat']));
    % parent grians of interest
    ID_p_iB_to_1_cell = d.ID_p_iB_to_1_cell;
    gList = d.gList;
    % some mask
    big_twin_grain_in_overlap_area_cell = d.big_twin_grain_in_overlap_area_cell;
    ID_overlap = d.ID_overlap;
    frd_twin_cell = d.frd_twin_cell;
    type_1_grain_ID_cell = d.type_1_grain_ID_cell;
    type_1_grain_label_cell = d.type_1_grain_label_cell;
    type_2_grain_ID_cell = d.type_2_grain_ID_cell;
    type_2_grain_label_cell = d.type_2_grain_label_cell;
    variant_pixel_mode = d.variant_pixel_mode;
    
    variable_names = {'iE','ID','iTwin','nvpxl_current','nvpxl_max','SF', 'was_active','detwinned'};
    tbl = cell2table(cell(0,length(variable_names)));
    tbl.Properties.VariableNames = variable_names;
    
    for iE = 4:7
        iB = iE + 1;
        % for parent grains, just use ID_overlap
        ID_p = ID_overlap;
        
        frw_twin = frd_twin_cell{iB};
        type_1_grain_ID = type_1_grain_ID_cell{iB};
        type_1_grain_label = type_1_grain_label_cell{iB};
        type_2_grain_ID = type_2_grain_ID_cell{iB};
        type_2_grain_label = type_2_grain_label_cell{iB};
        % mode_based vp_map (for type-1 grains): use variant_pixel_mode + mask_of_current_twin
        vp_map = variant_pixel_mode;
        vp_map(type_1_grain_ID==0) = 0;
        % mode_based vp_map for type-2 grains
        vp_map2 = variant_pixel_mode;
        vp_map2(~(type_2_grain_ID>0)) = 0;
        
        % [task-1 and 2]
        
        for ig = 1:length(gList)
            ID_current = gList(ig);
            
            % calculate variant SFs
            ind = find(gID==ID_current);
            euler_po = [gPhi1(ind), gPhi(ind), gPhi2(ind)];
            [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler_po, [0 0 0], [0 0 0], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
            tSFs = abs_schmid_factor(19:24,2);
            
            % Check if was active and if significantly detwinned:
            iv_was_active = zeros(6,1);
            iv_sig_detwinned = zeros(6,1);
            
            ind = ismember(ID_p, ID_current);
            t2_grains = type_2_grain_ID(ind);
            t2_grains = unique(t2_grains);
            t2_grains(t2_grains==0) = [];
            for it2 = 1:length(t2_grains)
                id_t2 = t2_grains(it2);
                ind = ismember(type_2_grain_ID, id_t2);
                
                % determine active variants at max compression
                v_pixels = vp_map2(ind);
                t2_variants = [];
                t2_label = unique(type_2_grain_label(ind)); % g4 or g5 grain ?     
                if numel(t2_label) > 1
                    error(' ');
                end
                for itwin = 1:6
                    if sum(v_pixels==itwin)/numel(v_pixels) > 0.1
                        t2_variants = [t2_variants, itwin];
                        if t2_label == 5
                            iv_sig_detwinned(itwin) = 1;
                        end
                    end
                end
                iv_was_active(t2_variants) = 1;
                
            end
            
            for iv = 1:6
                % current twin pixel
                ind = ismember(ID_overlap, ID_current) & ismember(vp_map, iv);
                nvpxl_current = sum(ind(:));
                % max twin pixel
                ind = ismember(ID_overlap, ID_current) & ismember(vp_map2, iv);
                nvpxl_max = sum(ind(:));
                % variant SF
                vsf = tSFs(iv);
                
                tbl = [tbl; {iE, ID_current, iv, nvpxl_current, nvpxl_max, vsf, iv_was_active(iv), iv_sig_detwinned(iv)}];
            end
        end
        
    end
    
    save(fullfile(output_dir, 'detwin', [sample_name,'_detwin.mat']), 'tbl');
end
close all;

%% [fig 9] [explore 1] boxplot: detwin fraction vs SF, each material, at iEs

edges = -0.5:0.05:0.5;
nGroups = length(edges)-1;
clear labels;
for ii = 1:length(edges)-1
    labels{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
end

for icell = [2,5,7,10] %1:size(cells,1)
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    sample_ID = cells{icell,3};
    sample_material = cells{icell,7};
    sample_label = cells{icell, 9};

    load(fullfile(output_dir, 'detwin', [sample_name,'_detwin.mat']), 'tbl');
    
    for iE = 4:7
        ind = (tbl.iE==iE)&(tbl.nvpxl_max>0); % twins
        vy = (tbl.nvpxl_max(ind) - tbl.nvpxl_current(ind))./tbl.nvpxl_max(ind);  % y = percent decrease
        vx  = tbl.SF(ind);  % vx = SFs
        gv = discretize(vx, edges);
        
        figure;disableDefaultInteractivity(gca);
        % boxplot([vy; nan*ones(nGroups, 1)], [gv; (1:nGroups)'],'notch','off');
        boxchart([gv; (1:nGroups)'], [vy; nan*ones(nGroups, 1)], 'linewidth',1.5);
        xlabel('Twin Variant Schmid Factor');
        ylabel('Variant Detwin Fraction')
        set(gca,'xTick',1:nGroups, 'xticklabels',labels,'xticklabelrotation',45,'fontsize',15, 'xlim',[8.5 20.5], 'ylim', [-0.1 1.1]);
        title_str = [sample_label, ', load step ',num2str(iE)];
        title(title_str,'fontweight','normal');
        print(fullfile(output_dir, 'detwin', [sample_name,'_by_SF_iE_',num2str(iE),'.tiff']), '-dtiff');
        close all;
    end
end

%% [explore 2, maybe not very useful] boxplot: detwin fraction by number_of_variant_pixels, each material, at iEs
close all;
edges = [0:200:3000,inf];
nGroups = length(edges)-1;
clear labels;
for ii = 1:length(edges)-1
    labels{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
end

for icell = [2,5,7,10] %1:size(cells,1)
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    sample_ID = cells{icell,3};
    sample_material = cells{icell,7};
    
    load(fullfile(output_dir, 'detwin', [sample_name,'_detwin.mat']), 'tbl');
    
    for iE = 4:7
        ind = (tbl.iE==iE)&(tbl.nvpxl_max>0); % twins
        vy = (tbl.nvpxl_max(ind) - tbl.nvpxl_current(ind))./tbl.nvpxl_max(ind);  % y = percent decrease
        vx  = tbl.nvpxl_max(ind);  % vx = max pixels of variant
        gv = discretize(vx, edges);
        
        figure;disableDefaultInteractivity(gca);
        boxplot([vy; nan*ones(nGroups, 1)], [gv; (1:nGroups)'],'notch','off');
        xlabel('# of Variant Pixels in Grain');
        ylabel('Variant Detwin Fraction')
        set(gca,'xticklabels',labels,'xticklabelrotation',45,'fontsize',15, 'ylim', [-0.1 1.1]);
        title_str = [sample_ID, ' iE=',num2str(iE)];
        title(title_str,'fontweight','normal');
        %         print(fullfile(output_dir, 'detwin', [sample_name,'_by_vol_iE_',num2str(iE),'.tiff']), '-dtiff');
        
    end
    close all;
end

%% [explore 3, maybe not very useful] This is silightly different than my previous 'significantly detwinned variant'.
edges = -0.5:0.05:0.5;
d_int = (edges(3)-edges(2))/2;
xpos = [edges(2)-d_int, edges(2:end-1)+d_int];

cs = inferno(5);
colors = [cs(1:4,:);
    cs(1:4,:);
    cs(1:3,:);
    cs(1:3,:)];

markers = {'-.o','-.d','-.s','-.^', ...
    '-.o','-.d','-.s','-.^', ...
    '-.o','-.d','-.s', ...
    '-.o','-.d','-.s', ...
    '-o','-d','-s','-^'};

close all;
for icell = [2,5,7,10] %1:size(cells,1)
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    sample_ID = cells{icell,3};
    sample_material = cells{icell,7};
    
    load(fullfile(output_dir, 'detwin', [sample_name,'_detwin.mat']), 'tbl');
    for iE = 4:7
        iB = iE + 1;
        ind = (tbl.iE==iE)&(tbl.was_active==1)&(tbl.detwinned==1);  % detwinned
        ind2 = (tbl.iE==iE)&(tbl.was_active==1)&(tbl.detwinned==0);  % twinned, not twinned
        
        [N_dt,~] = histcounts(tbl.SF(ind), edges);
        [N_ndt,~] = histcounts(tbl.SF(ind2), edges);

        
        figure; hold on; disableDefaultInteractivity(gca);
        hbar = bar(xpos, [N_ndt(:), N_dt(:)], 1, 'stacked');
        
        ymax = cells{icell,8};
        set(gca,'fontsize',12, 'XTick',-0.5:0.1:0.5, 'ylim', [0 ymax]);
        
        xlabel('Twin Variant Schmid Factor');
        ylabel('Counts');
        
        yyaxis right;
        set(gca, 'ycolor', 'k','fontsize',16, 'ylim',[0 150]);
        ylabel('Percent (%)');
        plot(-0.475:0.05:0.475, N_dt./(N_dt+N_ndt) * 100,'-ko','linewidth',1.5);
        pct(iB,:) = N_dt./(N_dt+N_ndt);
        
        title([sample_ID,', iE = ',num2str(iE)],'fontweight','normal');
        legend({'Variants Not Significantly Twinned', 'Variants Significantly Twinned','Pct of Variants Significantly Twinned'},'Location','northwest');
                
        print(fullfile(output_dir, 'detwin', [sample_name,'_sig_detwin_by_SF_iE_',num2str(iE),'.tiff']), '-dtiff');
    end
    
    % for all iEs
    figure; hold on;
    yyaxis right;
    hbar = bar(xpos, [N_ndt(:)+N_dt(:)], 1);
    set(gca,'ycolor', [.4 .4 .4]);
    set(gca, 'fontsize',13, 'XTick', -0.5:0.1:0.5, 'xlim',[-0.3 0.5]);
    hbar(1).FaceColor = [.85 .85 .85];
    xlabel('Twin Variant Schmid Factor');
    ylabel('Counts');
    
    yyaxis left;
    set(gca, 'ycolor', 'k', 'ylim',[0 100]);
    ylabel('Percent (%)');
    legend_str = [];
    for iE = 4:7
        iB = iE + 1;
        plot(edges(1:end-1)+0.025, 100*pct(iB,:), markers{iB}, 'color', colors(iB,:), 'linewidth', 1.5);
        legend_str = [legend_str, {['Load step ',num2str(iE)]}];
    end
    set(gca, 'ylim',[0 150], 'SortMethod', 'depth');
    legend(legend_str, 'Location','northwest');
    % print(fullfile(output_dir, 'detwin', [sample_name, ' sig detwin.tiff']),'-dtiff');

    close all;
end




