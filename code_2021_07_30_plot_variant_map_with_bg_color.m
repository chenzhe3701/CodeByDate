% make twin variant map, with good background.

working_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis';

output_dir = fullfile('E:\zhec umich Drive\0_temp_output\for 2021 workshop', 'variant maps');
mkdir(output_dir);

%% generate variant map
for iE = 0:13
    iB = iE + 1;
    % EBSD data for boundary file
    d = load(fullfile(working_dir, ['Mg4Al_U2_parent_grain_file_iE_',num2str(iE),'.mat']));
    
    x = d.x;
    y = d.y;
    ID = d.ID;
    
    boundary = find_one_boundary_from_ID_matrix(ID);
    
    % variant data
    d = load(fullfile(working_dir, 'variant_maps.mat'));
    if iE==0
        variant_pixel_map = zeros(size(ID));
    else
        variant_pixel_map = d.variant_point_wise{iE};
    end
    
    [f,a,c] = myplot(x,y,variant_pixel_map, boundary);
    
    make_variant_map_background;
    
    print(fullfile(output_dir, ['variant_pixel_map_iE_',num2str(iE),'.tiff']), '-dtiff');
    
    close all;
end

%%
winopen(output_dir);