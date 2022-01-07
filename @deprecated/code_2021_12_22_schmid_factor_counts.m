% Analyze Schmid factor distribution for twins in all samples
% ==> this (used data) that does not eliminate pre-existing twins. There are updated codes, so this one seems no longer needed.   

addChenFunction;
input_dir = 'E:\zhec umich Drive\0_temp_output\2021-12-22 analyze grains and twins';    % load twin pct from table.
output_dir = 'E:\zhec umich Drive\0_temp_output\2021-12-22 twin schmid factor';
mkdir(output_dir);

% [data_dir, sample_name, sample_ID, plot_symbol, group_number, ymax]
cells = {'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1, 300;
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1, 300;
    'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2, 350;
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2, 350;
    'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3, 200;
    'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3, 200;
    'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'o', 4, 450;
    'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 's', 4, 450;
    'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 'x', 4, 450;
    'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'o', 5, 200;
    'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 's', 5, 200;
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 'x', 5, 200};

%%
close all;

for ii = 1:size(cells,1)
    
    sample_dir = cells{ii,1};
    sample_name = cells{ii,2};
    sample_ID = cells{ii,3};

    % load variant maps at all load steps
    d = matfile(fullfile(sample_dir, 'variant_maps.mat'));
    variant_map_cell = d.variant_point_wise;
    % twin data
    d = matfile(fullfile(sample_dir, 'twin_pct.mat'));
    tbl_1 = d.tbl;
    strain_ebsd = tbl_1.strain_ebsd;
    
    % table to summarize twin area, corrected twin area, at each iE
    variable_names = {'iE', 'total_twin_pixels', 'total_twin_pixels_corrected', 'total_pixels', 'twin_pct_corrected', 'twin_strain_corrected'};
    tbl_3 = cell2table(cell(0,length(variable_names)));
    tbl_3.Properties.VariableNames = variable_names;    
    
    % reference parent orientation data from iE=0
    d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_0.mat']));
    ID = d.ID;
    gPhi1 = d.gPhi1;
    gPhi = d.gPhi;
    gPhi2 = d.gPhi2;
    gID = d.gID;
    gDiameter = d.gDiameter;
    gArea = d.gArea;
    
    % analyze each load step iE
    for iE = 0 : length(variant_map_cell)-1
        iB = iE + 1;
        
        variant_map = variant_map_cell{iB};
        
        % load tbl_2, where  variable_names = {'iE','ID','iTwin','tSF', 'n_variant_pixels', 'variant_area_corrected', 'gDiameter','gArea'};
        load(fullfile(input_dir, [sample_name,'_iE=',num2str(iE),'.mat']),'tbl_2');

        % parent grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        % boundary_p = find_one_boundary_from_ID_matrix(ID_p);
        
        % child grain data
        % d = matfile(fullfile(sample_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));
               
        edges = -0.5:0.05:0.5;

        ind = (tbl_2.iE==iE)&(tbl_2.n_variant_pixels>0);    % twinned variants
        ind2 = (tbl_2.iE==iE)&(tbl_2.n_variant_pixels==0);
        [N_t,~] = histcounts(tbl_2.tSF(ind), edges);
        [N_nt,~] = histcounts(tbl_2.tSF(ind2), edges);
        
        figure; hold on; disableDefaultInteractivity(gca);
        hbar = bar(edges(1:end-1)+0.025, [N_nt(:), N_t(:)], 1, 'stacked');
        
        ymax = cells{ii,6};
        set(gca,'fontsize',12, 'XTick',-0.5:0.1:0.5,'ylim',[0 ymax]);

        xlabel('Twin Variant Schmid Factor');
        ylabel('Counts');
        
        yyaxis right;
        set(gca, 'ycolor', 'k','fontsize',16, 'ylim',[0 100]);
        ylabel('Percent (%)');
        plot(-0.475:0.05:0.475, N_t./(N_t+N_nt) * 100,'-ko','linewidth',1.5);
        
        title([sample_ID,', iE = ',num2str(iE)],'fontweight','normal');
        legend({'Variants Not Twinned', 'Variants Twinned','Percent of Variants Twinned'},'Location','northwest');
        
        print(fullfile(output_dir, [sample_ID,'  iE = ',num2str(iE),'.tiff']),'-dtiff');
        close all;
    end    


end
