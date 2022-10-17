% check at each load step, # of each type-1 type-2 grains

clear; clc; close all;
addChenFunction;

% [data_dir, sample_name, sample_ID, plot_symbol, group_number, ymax, sample_material]
cells = {'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1, 50, 'Mg4Al FG', 'Mg4Al fine grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1, 50, 'Mg4Al FG', 'Mg4Al fine grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2, 350, 'Mg4Al MG', 'Mg4Al medium grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2, 350, 'Mg4Al MG', 'Mg4Al medium grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3, 100, 'Mg4Al CG', 'Mg4Al coarse grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3, 100, 'Mg4Al CG', 'Mg4Al coarse grain';
    % 'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'x', 4, 30, 'Mg FG', 'Mg fine grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 'o', 4, 30, 'Mg FG', 'Mg fine grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 's', 4, 30, 'Mg FG', 'Mg fine grain';
    % 'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'x', 5, 200, 'Mg CG', 'Mg coarse grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 'o', 5, 200, 'Mg CG', 'Mg coarse grain';
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 's', 5, 200, 'Mg CG', 'Mg coarse grain'};

% variant_map_dir = 'H:\Other computers\My Laptop w541 2022\zhec umich Drive\All twin variant maps cleaned';

% location of the twin evolution (twin detwin retwin) data
evolution_dir = 'C:\Users\chenz\Work\Data\0_temp_output\all twin evolution analysis';

output_dir = 'C:\Users\chenz\Work\Data\0_temp_output\all twin evolution analysis';
mkdir(output_dir);

%%
mkdir(fullfile(output_dir, 'type gn evolution'));     
for icell = [2,5,7,10] % 1:size(cells,1)
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    sample_ID = cells{icell,3};
    sample_material = cells{icell,7};
    sample_label = cells{icell, 8};
        
    % twin evolution data
    d = matfile(fullfile(evolution_dir, [sample_name, ' twin evolution label.mat']));
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
    
    myMat = zeros(14,5);
    for iE = 0:13
        iB = iE + 1;
        % for parent grains, just use ID_overlap
        ID_p = ID_overlap;
        
        %         frw_twin = frd_twin_cell{iE};
        type_1_grain_ID = type_1_grain_ID_cell{iB};
        type_1_grain_label = type_1_grain_label_cell{iB};
        type_2_grain_ID = type_2_grain_ID_cell{iB};
        type_2_grain_label = type_2_grain_label_cell{iB};
        
        g1_new = 1; % (g1) new twin type-1 grain
        g2_detwin_retwin = 2; % (g2) detwin then retwin type-1 grain
        g3_evolving_1 = 3; % (g3) evolving type-1 grain
        % type-2 grain, past-or-present twin grain (ever twinned up to iE)
        g4_evolving_2 = 4; % (g4) evolving type-2 twin
        g5_detwin_2 = 5; % (g5) completely detwin type-2 grain
        
        for iLabel = 1:3
            ind = ismember(type_1_grain_label, iLabel);
            ngs = unique(type_1_grain_ID(ind));
            if any(ngs==0)
                error(' ');
            end
            if ~isempty(ngs)
                myMat(iB,iLabel) = numel(ngs);
            end
        end
        for iLabel = 4:5
            ind = ismember(type_2_grain_label, iLabel);
            ngs = unique(type_2_grain_ID(ind));
            if any(ngs==0)
                error(' ');
            end
            if ~isempty(ngs)
                myMat(iB,iLabel) = numel(ngs);
            end
        end
        
    end
    
    colors = plasma(25);
    g_color{1} = colors(24,:);    % g1 = completely new twin, yellowish
    g_color{2} = colors(10,:);     % g2 = de-twin then re-twin, purple
    g_color{3} = colors(16,:);     % g3 = evoling current twin, pink
    g_color{4} = [0.7, 0.7, 0.7];     % g4 = evolving past-or-present twin, light gray
    g_color{5} = [0.2, 0.2, 0.2];     % g5 = completely detwinned twin, blk 
    
    g_str = {'g1: New Twin Type-1 Grain', ...
        'g2: Detwin Then Retwin Type-1 Grain', ...
        'g3: Evolving Type-1 Grain', ...
        'g4: Evolving Type-2 Grain', ...
        'g5: Completely Detwin Type-2 Grain'};
    
    close all; 
    figure;
    set(gcf,'Position', [100 100 900, 450]);
    types_to_plot = [3,2,1,5];
    hbar = bar(0:13, myMat(:,types_to_plot), 0.9, 'stacked');
    for kk = 1:length(types_to_plot)
        hbar(kk).FaceColor = g_color{types_to_plot(kk)};
    end
    set(gca,'fontsize',16, 'XTickLabelRotation', 0);
    title(sample_label, 'fontweight','normal');
    xlabel('Load Step');
    ylabel('Number of Grains');
    legend(g_str(types_to_plot), 'Location', 'eastoutside');
    
    print(fullfile(output_dir, 'type gn evolution', [sample_name,' type gn evltn.tiff']), '-dtiff');
end











