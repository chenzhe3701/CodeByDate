% combine multiple images into a single image for each load step

output_dir = 'E:\zhec umich Drive\2021-08-20 UM129 Mg_C2 insitu EBSD\analysis\merged';
mkdir(output_dir);
%%  Dir for resources. Template for file names.

% sources of variant map, experimentally measured position
vMap_dir = 'E:\zhec umich Drive\2021-08-20 UM129 Mg_C2 insitu EBSD\analysis\variant maps';
vMap_file = ['variant_pixel_map_iE_',num2str(iE),'.tiff'];

% spatially deformed to iE=0 variant map
vMap_spatially_dir = 'E:\zhec umich Drive\0_temp_output\2021-05-04 Analyze persistent twin\Mg4Al_U2';
vMap_spatially_file = ['pixel variant map iE=',num2str(iE),'.tiff'];

IPF_dir = 'E:\zhec umich Drive\2021-08-20 UM129 Mg_C2 insitu EBSD';
IPF_file = ['UM129_Mg_C2 parent IPF ND iE=',num2str(iE),'.tif'];

% cat_dir = 'E:\zhec umich Drive\0_temp_output\2021-05-04 Analyze persistent twin\Mg4Al_U2';
% cat_file = ['frd twinned iE=',num2str(iE),'.tiff'];     % fresh, re, de-twins 

curve_dir = 'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu curve';
curve_file = ['stress vs strain iE_',num2str(iE),'.tiff'];

%% For each iE, merge the maps
for iE = 0:13
    close all;
    % (1) Loading curve
    curve_file = ['stress vs strain iE_',num2str(iE),'.tiff'];
    image_1 = imread(fullfile(curve_dir, curve_file));
    % (2) IPF
    IPF_file = ['UM129_Mg_C2 parent IPF ND iE=',num2str(iE),'.tif'];
    image_2 = imread(fullfile(IPF_dir, IPF_file));
    % (3) twin variant map
    vMap_file = ['variant_pixel_map_iE_',num2str(iE),'.tiff'];
    image_3 = imread(fullfile(vMap_dir, vMap_file));
    % (4) category map
%     cat_file = ['frd twinned iE=',num2str(iE),'.tiff']; 
%     image_4 = imread(fullfile(cat_dir, cat_file));
    
%     figure; imshow(image_1);
%     figure; imshow(image_2);
%     figure; imshow(image_3);
%     figure; imshow(image_4);

    I_1 = imcrop(image_1, [0,0, inf,inf]);
    I_1 = imresize(I_1, 620/656);
    
    I_2 = imcrop(image_2, [0,0, inf,inf]);
    I_2 = imresize(I_2, [761, 761]);
    
    I_3 = imcrop(image_3, [200,70, 880,765]);
    I_3 = imresize(I_3, 761/765);
    
%     I_4 = imcrop(image_4, [215,75, inf,760]);
    
    I = inf * uint8(ones(1500,1900,3));
    
    % stitch I_1
    indR = 30;
    indC = 950;
    I(indR:indR + size(I_1,1) - 1, indC:indC + size(I_1,2) - 1, :) = I_1;
    
    % stitch I_2
    indR = 710;
    indC = 50;
    I(indR:indR + size(I_2,1) - 1, indC:indC + size(I_2,2) - 1, :) = I_2;
    
    % stitch I_3
    indR = 710;
    indC = 950;
    I(indR:indR + size(I_3,1) - 1, indC:indC + size(I_3,2) - 1, :) = I_3;
    
    % stitch I_4
%     indR = 710;
%     indC = 1900;
%     I(indR:indR + size(I_4,1) - 1, indC:indC + size(I_4,2) - 1, :) = I_4;
    
    figure;
    imshow(I);
    
    % print(fullfile(output_dir, ['print merged_iE=',num2str(iE),'.tiff']),'-dtiff');    
    imwrite(I, fullfile(output_dir, ['merged_iE=',num2str(iE),'.tiff']));
end

%%
close all;



