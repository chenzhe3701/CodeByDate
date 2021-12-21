% check if 'each grain in the child grain data cannot belong to more than
% one grains in the parent grain data'

cells = {'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'o', 1;
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 's', 1;
    'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'o', 2;
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 's', 2;
    'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'o', 3;
    'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 's', 3;
    'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'o', 4;
    'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 's', 4;
    'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'x', 4;
    'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'o', 5;
    'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 's', 5;
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'x', 5};
%%

for ii = 1:length(cells)
    sample_dir = cells{ii,1};
    sample_name = cells{ii,2};
    
    for iE = 0:13
        iB = iE + 1;
        
        disp(sample_name);
        disp(iE);
        
        % parent grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        
        % child grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));       
        ID_c = d.ID;
        
        uniqueIDc = unique(ID_c(:));
        for jj = 1:length(uniqueIDc)
            IDc_current = uniqueIDc(jj);
            inds = ismember(ID_c, IDc_current);
            
            parent_IDs = unique(ID_p(inds));
            if length(parent_IDs) > 1
                error(' ');
            end
        end
    end
end