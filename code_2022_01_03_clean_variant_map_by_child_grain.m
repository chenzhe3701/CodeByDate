% clean up twin maps
% (1) Rmove twin identified for child grains < 16 data pts  

addChenFunction;
output_dir = 'E:\zhec umich Drive\All twin variant maps cleaned';
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
    
    % load twin maps
    d = matfile(fullfile(sample_dir,'variant_maps.mat'));
    variant_point_wise = d.variant_point_wise;
    variant_grain_wise = d.variant_grain_wise;
    
    % run this for ID_list{iB}
    run(fullfile(correction_dir, ['variables_',sample_name,'.m']));
    
    for iE = 0:13
        iB = iE + 1;
        
        % load unmodified parent grain data, find parent grains that were modified
        d = matfile(fullfile(sample_dir, 'step-1', [sample_name, '_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        
        % mask for modified parent grain vs. unmodified parent grain
        mask_mod = ismember(ID_p, ID_list{iB});
        bd = find_one_boundary_from_ID_matrix(mask_mod+1);
        
        vp = variant_point_wise{iB};
        vg = variant_grain_wise{iB};
        myplot(vp, bd);
        print(fullfile(output_dir, [sample_name,'_vp_iE_',num2str(iE),'_old.tif']),'-dtiff');
        myplot(vg);
        print(fullfile(output_dir, [sample_name,'_vg_iE_',num2str(iE),'_old.tif']),'-dtiff');
        
        % load child grain
        d = matfile(fullfile(sample_dir, [sample_name, '_grain_file_iE_',num2str(iE),'.mat']));
        ID_c = d.ID;
        id_c_list = []; % id_c that contains small new grains
        % find unique child grain list
        % ID_c_list = unique(ID_c(:));
        
        % one_pass_label to check child grain size
        ID_c_temp = one_pass_label(ID_c);
        ID_c_temp_list = unique(ID_c_temp(:));

        for jj = 1:length(ID_c_temp_list)
            id_c_temp = ID_c_temp_list(jj);
            disp([sample_name,' iE=',num2str(iE),',iC=',num2str(id_c_temp)]);
            
            inds = ismember(ID_c_temp, id_c_temp);
            
            gs = sum(inds(:));
            if gs<16
                vp(inds) = 0;
                vg(inds) = 0;
                
                id_c = ID_c(inds);  % this area belong to the original id_c   
                id_c_list = unique([id_c_list; id_c]);
            end
        end
        
        % copy back
        variant_point_wise{iB} = vp;
       	variant_grain_wise{iB} = vg;
        
        myplot(vp);
        print(fullfile(output_dir, [sample_name,'_vp_iE_',num2str(iE),'.tif']),'-dtiff');
        myplot(vg);
        print(fullfile(output_dir, [sample_name,'_vg_iE_',num2str(iE),'.tif']),'-dtiff');
        
        ID_c(~ismember(ID_c,id_c_list)) = 0;
        myplot(ID_c); 
        try
            caxis([min(id_c_list)-1, max(id_c_list)+1]);
        end
        print(fullfile(output_dir, [sample_name,'_IDcAffected_iE_',num2str(iE),'.tif']),'-dtiff');
        close all;
    end
    
    save(fullfile(output_dir, [sample_name, '_variant_maps.mat']), 'variant_point_wise', 'variant_grain_wise');
end




