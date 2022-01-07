% I would like to check the 'child grain map', and plot grains with small
% number of pixels.  They are likely introduced in the grain modification
% process.
% Answer: yes. Acutally, depending on the clean up process, anything < 16
% data points should be artifact, and should not be considered as twin.

addChenFunction;
output_dir = 'E:\zhec umich Drive\0_temp_output\2022-01-03 temp';
mkdir(output_dir);

% [data_dir, sample_name, sample_ID, plot_symbol, group_number]
cells = {'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1;
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1;
    'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2;
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2;
    'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3;
    'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3;
    'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'x', 4;
    'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 'o', 4;
    'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 's', 4;
    'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'x', 5;
    'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 'o', 5;
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 's', 5};

correction_dir = 'D:\p\m\Twin_Identification_from_EBSD\setting files for samples';

%%
close all;
for ii = 1:12
    
    sample_dir = cells{ii,1};
    sample_name = cells{ii,2};
    sample_ID = cells{ii,3};
    
    % run this for ID_list{iB}
    run(fullfile(correction_dir, ['variables_',sample_name,'.m']));
    
    % load sample setting file, look for which (original) parent grain is modified
    for iE = 0:13
        iB = iE + 1;
        % load unmodified parent grain data, find parent grains that were modified
        d = matfile(fullfile(sample_dir, 'step-1', [sample_name, '_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
                
        % mask for modified parent grain vs. unmodified parent grain
        mask_mod = ismember(ID_p, ID_list{iB});
        bd = find_one_boundary_from_ID_matrix(mask_mod+1);
        
        % load child grain data
        d = matfile(fullfile(sample_dir, [sample_name, '_grain_file_iE_',num2str(iE),'.mat']));
        ID_c = d.ID;

        % for each child grain, count grain size
        ID_c_list = unique(ID_c(:));
        
        % Initialize
        child_gs_map = zeros(size(ID_c));
        id_size_under_mask = [];
        id_size_out_mask = [];
        
        for jj = 1:length(ID_c_list)
           id_c = ID_c_list(jj); 
           inds = ismember(ID_c, id_c);

           gs = sum(inds(:));
           child_gs_map(inds) = gs;
           if sum(mask_mod(inds))>0
               % child grain in modified area
               id_size_under_mask = [id_size_under_mask; id_c, gs]; 
           else
               id_size_out_mask = [id_size_out_mask; id_c, gs]; 
           end
           
                     
        end
        
        figure; 
        histogram(id_size_under_mask(:,2), 1:100);
        print(fullfile(output_dir, [sample_name, '_cgs_mod_iE_',num2str(iE),'.tif']), '-dtiff');
        figure;
        histogram(id_size_out_mask(:,2), 1:100);
        set(gca,'XMinorTick','on')
        print(fullfile(output_dir, [sample_name, '_cgs_nomod_iE_',num2str(iE),'.tif']), '-dtiff');
        
        myplot(child_gs_map<16, bd);
        print(fullfile(output_dir, [sample_name, '_small_cGrain_mod_iE_',num2str(iE),'.tif']), '-dtiff');
        
        close all;
    end
    
end







