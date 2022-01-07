%% this code is replaced by the paper fig SF distribution code.
% Study the twin contributed strain (corrected by twin SF) vs.
% (o) global strain estimated from EBSD
% (o) twin area fraction

% Method for calculate twin contributed strain:
% twin_contributed_strain = 0.129 * sum_over_grains *
% sum_over_variants(num_of_variant_pixels * SF) / total_num_of_twin_pixels

% Note that the processed data follows the OIM convention for grain
% orientation.  For a grain mostly twinned, the grain orientation is the
% twin orientation.  Therefore, we always need to load data of iE=0 as a
% reference for the parent orientation!!!

addChenFunction;
input_dir = 'E:\zhec umich Drive\0_temp_output\2021-12-22 analyze grains and twins';
output_dir = 'E:\zhec umich Drive\0_temp_output\2021-12-22 negative SF twins';
mkdir(output_dir);

% [data_dir, sample_name, sample_ID, plot_symbol, group_number]
cells = {'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1;
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1;
    'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2;
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2;
    'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3;
    'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3;
    'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'o', 4;
    'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 's', 4;
    'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 'x', 4;
    'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'o', 5;
    'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 's', 5;
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 'x', 5};

%%
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
        
    % reference parent orientation data from iE=0
    d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_0.mat']));
    ID = d.ID;
    gPhi1 = d.gPhi1;
    gPhi = d.gPhi;
    gPhi2 = d.gPhi2;
    gID = d.gID;
    gDiameter = d.gDiameter;
    gArea = d.gArea;
    
    uniqueID = [];    
    % First need to go through all iEs, find unique parent grains to summarize
    for iE = 0:13
        % parent grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        gID_p = d.gID;
        if iE == 0 
            uniqueID = gID_p;
        else
            uniqueID = intersect(uniqueID, gID_p);
        end
    end
    
    % analyze each load step iE
    for iE = 0:13
        iB = iE + 1;
        
        variant_map = variant_map_cell{iB};
        
        % with variable_names = {'iE','ID','iTwin','tSF', 'n_variant_pixels', 'variant_area_corrected', 'gDiameter','gArea'};
        load(fullfile(input_dir, [sample_name,'_iE=',num2str(iE),'.mat']),'tbl_2');

        % parent grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        boundary_p = find_one_boundary_from_ID_matrix(ID_p);
        
        % child grain data
        % d = matfile(fullfile(sample_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));
        
        for jj = 1:length(uniqueID)
            ID_current = uniqueID(jj);
            
            for iTwin = 1:6
                indt = (tbl_2.ID==ID_current)&(tbl_2.iTwin==iTwin);
                tSF = tbl_2.tSF(indt);
                
                % if positive twin SF, remove from the 'variant map'
                if tSF > 0
                    ind_twin = ismember(ID_p, ID_current) & ismember(variant_map, iTwin);
                    variant_map(ind_twin) = 0;
                end
            end            
        end
        
        % after checking all grains, plot and save
        myplot(variant_map, boundary_p);
        print(fullfile(output_dir, [sample_ID,'  iE = ',num2str(iE),'.tiff']),'-dtiff');
        close all;
    end

end



