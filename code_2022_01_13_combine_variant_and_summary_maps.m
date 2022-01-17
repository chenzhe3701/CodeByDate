% make evolution maps, combining images
% (1) variant map, (2) pixel summary, (3) grain summary

clear; clc; close all;
addChenFunction;

% [data_dir, sample_name, sample_ID, plot_symbol, group_number, ymax]
cells = {'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1, 300;
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1, 300;
    'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2, 350;
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2, 350;
    'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3, 200;
    'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3, 200;
    % 'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'x', 4, 450;
    'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 'o', 4, 450;
    'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 's', 4, 450;
    % 'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'x', 5, 200;
    'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 'o', 5, 200;
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 's', 5, 200};

variant_map_dir = 'E:\zhec umich Drive\All twin variant maps cleaned';
variant_img_dir = 'E:\zhec umich Drive\0_temp_output\variant maps bgbk\crop to only map';

evolution_img_dir = 'E:\zhec umich Drive\0_temp_output\all twin evolution analysis\evolution map';

output_dir = 'E:\zhec umich Drive\0_temp_output\all twin evolution analysis\combined maps';
mkdir(output_dir);
%%

for icell = 1:10
    
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    sample_ID = cells{icell,3};
    
    for iE = 0:13
        iB = iE + 1;
        % (1) Variant map
        image_1 = imread(fullfile(variant_img_dir,  sample_name, [sample_name, 'variant_pixel_map_iE_',num2str(iE),'.tiff']));
        % (2) pixel level summary map
        image_2 = imread(fullfile(evolution_img_dir, [sample_name, ' frd twinned iE=',num2str(iE),'.tiff']));
        % (3) grain level summary map
        image_3 = imread(fullfile(evolution_img_dir, [sample_name, ' grain label iE=',num2str(iE),'.tiff']));
        
        % figure; imcrop(image_3);
        
        I_1 = imcrop(image_1, [1 1 766 766]);
        I_1 = imresize(I_1, 750/766);
        
        I_2 = imcrop(image_2, [379 70 766 766]);
        I_2 = imresize(I_2, 750/766);
        I_2([1,end],:,:) = 1;
        I_2(:,[1,end],:) = 1;
        
        I_3 = imcrop(image_3, [541 70 766 766]);
        I_3 = imresize(I_3, 750/766);
        I_3([1,end],:,:) = 1;
        I_3(:,[1,end],:) = 1;
        
        I = inf * uint8(ones(770,2370,3));
        % stitch I_1
        indR = 10;
        indC = 10;
        I(indR:indR + size(I_1,1) - 1, indC:indC + size(I_1,2) - 1, :) = I_1;
        
        % stitch I_2
        indC = 810;
        I(indR:indR + size(I_2,1) - 1, indC:indC + size(I_2,2) - 1, :) = I_2;
        
        % stitch I_3
        indC = 1610;
        I(indR:indR + size(I_3,1) - 1, indC:indC + size(I_3,2) - 1, :) = I_3;
        
        % stitch I_4
        %     indR = 710;
        %     indC = 1900;
        %     I(indR:indR + size(I_4,1) - 1, indC:indC + size(I_4,2) - 1, :) = I_4;
        
        figure;
        imshow(I);
        
        imwrite(I, fullfile(output_dir, [sample_name,' iE=',num2str(iE),'.tiff']));
        close all;
    end
end

