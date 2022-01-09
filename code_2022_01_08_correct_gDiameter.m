% after data processing, I realized that the gDiameter calculated in and after
% step-2 was wrong.  They were incorrect as sqrt(4*gArea).  They should be
% sqrt(4 * gArea / pi).
clear; clc;

cells = {'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1, 300;
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1, 300;
    'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2, 350;
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2, 350;
    'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3, 200;
    'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3, 200;
    'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'x', 4, 450;
    'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 'o', 4, 450;
    'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 's', 4, 450;
    'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'x', 5, 200;
    'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 'o', 5, 200;
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 's', 5, 200;
    'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD\analysis', 'Mg4Al_C1', 'Mg4Al C1', 'x', 1, 300};

for icell = 1:size(cells,1)
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    for iE = 0:13
        % local dirs
        local_dirs = {'','step-1', 'step-2','step-3','step-5'};
        
        for jj = 1:5
            local_dir = local_dirs{jj};
            clear gDiameter gArea
            try
                fname = fullfile(sample_dir, local_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']);
                disp(fname);
                load(fname, 'gDiameter', 'gArea');
                if (pi/4*(gDiameter(1)).^2 - gArea(1)) < 0.001
                    disp(['OK: ', fname]);
                else
                    % this should be corrected
                    if 1/4*(gDiameter(1)).^2 - gArea(1) > 0.001
                        disp(['error: icell=', num2str(icell), 'iE=',num2str(iE)]);
                    end
                    gDiameter = (gArea*4/pi).^(1/2);
                end
                save(fname, 'gDiameter','-append');
            end
            
            clear gDiameter gArea
            try
                fname = fullfile(sample_dir, local_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']);
                disp(fname);
                load(fname, 'gDiameter', 'gArea');
                if (pi/4*(gDiameter(1)).^2 - gArea(1)) < 0.001
                    disp(['OK: ', fname]);
                else
                    % this should be corrected
                    if 1/4*(gDiameter(1)).^2 - gArea(1) > 0.001
                        disp(['error: icell=', num2str(icell), 'iE=',num2str(iE)]);
                    end
                    gDiameter = (gArea*4/pi).^(1/2);
                end
                save(fname, 'gDiameter','-append');
            end
            
        end
    end
end