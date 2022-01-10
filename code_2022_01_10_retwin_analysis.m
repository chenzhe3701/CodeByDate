% analyze detwin

clear; clc; close all;
addChenFunction;

% [data_dir, sample_name, sample_ID, plot_symbol, group_number, ymax, sample_material]
cells = {'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1, 50, 'Mg4Al FG';
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1, 50, 'Mg4Al FG';
    'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2, 350, 'Mg4Al MG';
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2, 350, 'Mg4Al MG';
    'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3, 100, 'Mg4Al CG';
    'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3, 100, 'Mg4Al CG';
    % 'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'x', 4, 30, 'Mg FG';
    'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 'o', 4, 30, 'Mg FG';
    'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 's', 4, 30, 'Mg FG';
    % 'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'x', 5, 200, 'Mg CG';
    'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 'o', 5, 200, 'Mg CG';
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 's', 5, 200, 'Mg CG'};

% variant_map_dir = 'E:\zhec umich Drive\All twin variant maps cleaned';

% location of the twin evolution (twin detwin retwin) data
evolution_dir = 'E:\zhec umich Drive\0_temp_output\all twin evolution analysis';

output_dir = 'E:\zhec umich Drive\0_temp_output\all twin evolution analysis';
mkdir(output_dir);

%% [analyze] retwin
mkdir(fullfile(output_dir, 'retwin'));

for icell = 1:size(cells,1)
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    sample_ID = cells{icell,3};
    sample_material = cells{icell,7};
    
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
    
    variable_names = {'iE','ID','iTwin','nvpxl_current','nvpxl_max','SF'};
    tbl = cell2table(cell(0,length(variable_names)));
    tbl.Properties.VariableNames = variable_names;
    
    N = [];
    for iE = 7:10
        iB = iE + 1;
        % for parent grains, just use ID_overlap
        ID_p = ID_overlap;
        
        frw_twin = frd_twin_cell{iE};
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
        
        % algorithm:
        %  frd_twin @ iE_7, divied into re-taining and detwin
        % For each load step, compare to load step 7
        %  accumulative new twin area in cycle 2 = (type-1 map) not in (type-2 map @ iE 7) 
        %  recurring twin area in cycle 2 = (type-1 map) in (type-2 map @ iE 7)   
        %  not-reactivated detwin area in cycle 2 = (type-2 map @ iE 7) not in (type-1 map)     
        
        % the reference step is ie=7
        ie = 7;
        ib = ie + 1;
        ind_1 = (type_1_grain_ID_cell{iB}>0) & ~(type_2_grain_ID_cell{ib}>0);
        ind_2 = (type_1_grain_ID_cell{iB}>0) & (type_2_grain_ID_cell{ib}>0);
        ind_3 = (type_2_grain_ID_cell{ib}>0) & ~(type_1_grain_ID_cell{iB}>0);
        N = [N; [sum(ind_3(:)), sum(ind_2(:)), sum(ind_1(:))]];
        
    end

    save(fullfile(output_dir, 'retwin', [sample_name,'_retwin.mat']), 'N');
end

%% [plot types of pixels in retwin process]
mkdir(fullfile(output_dir, 'retwin'));

for icell = [2,5,7,10] % 1:size(cells,1)
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    sample_ID = cells{icell,3};
    sample_material = cells{icell,7};
    
    load(fullfile(output_dir, 'retwin', [sample_name,'_retwin.mat']), 'N');
       
    figure;
    set(gcf,'Position', get(gcf,'Position') .* [1,1,0,0] + [0,0,850,420]);
    hbar = bar(N, 'stacked');
    colors = parula(16); 
    hbar(1).FaceColor = [.5 .5 .5];
    hbar(2).FaceColor = colors(2,:);
    hbar(3).FaceColor = colors(14,:);
    ymax = get(gca,'ylim');
    ymax = ymax(2);
    set(gca,'fontsize',14, 'XTickLabel', {'7','8','9','10'}, 'ylim', [0 ymax]);
    xlabel('Load Step');
    ylabel('Pixel Counts');
    
    yyaxis right;
    set(gca,'ycolor','k','ylim',[0, ymax/sum(N(1,:))]*100);
    ylabel('Percentage (%)');
    
    legend({'Detwin pixels', 'Recurring twin pixels', 'Accumulative new twin pixels'}, 'location','eastoutside');
    print(fullfile(output_dir, 'retwin', [sample_name,'_retwin.tiff']), '-dtiff');
    close all;
    
end

