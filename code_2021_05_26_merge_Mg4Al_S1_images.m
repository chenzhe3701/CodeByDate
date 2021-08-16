% For the twin-detwin-retwin maps, I want to put pixel level summary on the
% left, and grain level summary on the right
% Use this code to stitch the images

working_dir = 'E:\zhec umich Drive\0_temp_output\Mg4Al_S1 analysis';
save_dir = fullfile(working_dir, 'merged image');
mkdir(save_dir);

%% merge left and right
for iE = 2:6
    
    image_1 = imread(fullfile(working_dir, ['frd twinned iE=',num2str(iE),'.tiff']));
    image_2 = imread(fullfile(working_dir, ['grain label iE=',num2str(iE),'.tiff']));
    
    % [I, rectout] = imcrop(image_1);
    
    I_1 = imcrop(image_1, [100, 20, 1100, 850]);
    I_2 = imcrop(image_2, [290, 20, 1240, 850]);
    I = [I_1, I_2];
    
    imwrite(I, fullfile(save_dir, ['2 merged_iE=',num2str(iE),'.tiff']));
    
end

%% merge 4 images

for iE = 2:6
    close;
    image_1 = imread(fullfile(working_dir, ['exx_map_iE=',num2str(iE),'.tiff']));
    image_2 = imread(fullfile(working_dir, ['variant_map_iE=',num2str(iE),'.tiff']));
    image_3 = imread(fullfile(working_dir, ['frd twinned iE=',num2str(iE),'.tiff']));
    image_4 = imread(fullfile(working_dir, ['grain label iE=',num2str(iE),'.tiff']));
    
    % [I_1, rectout_a] = imcrop(image_1);
    
    I_1 = imcrop(image_1, [100, 20, 1100, 850]);
    I_2 = imcrop(image_2, [125, 20, 1100, 850]);
    I_3 = imcrop(image_3, [100, 20, 1100, 850]);
    I_4 = imcrop(image_4, [290, 20, 1240, 850]);
    I = [I_1, I_2, inf*ones(851,140,3); I_3, I_4];
    
    imwrite(I, fullfile(save_dir, ['4 merged_iE=',num2str(iE),'.tiff']));
    
end