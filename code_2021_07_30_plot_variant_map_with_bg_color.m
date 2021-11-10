%% make twin variant map, with good background, for all samples.

working_dirs = {% 'E:\Mg4Al_S1_insitu\EBSD Data\Mg4Al s1.osc', ...
    'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD\analysis', ...
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', ...
    'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', ...
    'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', ...
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', ...
    'E:\zhec umich Drive\2020-12-05 UM134 Mg_C1 insitu EBSD\analysis', ...
    'E:\zhec umich Drive\2021-01-15 UM134 Mg_C2 insitu EBSD\analysis', ...
    'E:\zhec umich Drive\2021-01-29 UM134 Mg_C3 insitu EBSD\analysis', ...
	'E:\zhec umich Drive\2021-06-29 UM129 Mg_C1 insitu EBSD\analysis', ...
    'E:\zhec umich Drive\2021-08-20 UM129 Mg_C2 insitu EBSD\analysis', ...
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis'};
sample_names = {'Mg4Al_C1', 'Mg4Al_C3', 'Mg4Al_U2', 'Mg4Al_A1', 'Mg4Al_A2', ...
    'UM134_Mg_C1', 'UM134_Mg_C2', 'UM134_Mg_C3', 'UM129_Mg_C1', 'UM129_Mg_C2', 'UM129_Mg_C3'};

output_dir = fullfile('E:\zhec umich Drive\0_temp_output\variant maps bgbk');

%% generate variant map
for ii = 1:length(sample_names)
    % make an output dir for this sample
    mkdir(fullfile(output_dir, sample_names{ii}));
    
    for iE = 0:13
        iB = iE + 1;
        
        % EBSD data for boundary file
        grain_file = fullfile(working_dirs{ii}, [sample_names{ii}, '_parent_grain_file_iE_',num2str(iE),'.mat']);
        if ~exist(grain_file,'file')
            break;
        end
        
        d = load(grain_file);
        
        x = d.x;
        y = d.y;
        ID = d.ID;
        
        boundary = find_one_boundary_from_ID_matrix(ID);
        
        % variant data
        d = load(fullfile(working_dirs{ii}, 'variant_maps.mat'));
        if iE==0
            variant_pixel_map = zeros(size(ID));
        else
            variant_pixel_map = d.variant_point_wise{iE};
        end
        
        [f,a,c] = myplot(x,y,variant_pixel_map, boundary);
        
        make_variant_map_background;
        
        print(fullfile(output_dir, sample_names{ii}, [sample_names{ii}, 'variant_pixel_map_iE_',num2str(iE),'.tiff']), '-dtiff');
        
        close all;
    end
end

%%
winopen(output_dir);

