working_dir = 'E:\zhec umich Drive\2021-08-30 try external in Dow';
hw = 100;
output_dir = fullfile(working_dir, 'analysis');
mkdir(output_dir);

%% [A] Increasing dwell time can lead to improved image quality
close all;

I = imread(fullfile(working_dir, 'ref_s4_r1_40um_d0.16.tiff'));
I = imcrop(I, [1885-hw, 2050-hw, 2*hw, 2*hw]);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, 'A_dwell_4us.tiff'));

I = imread(fullfile(working_dir, 'def_s8_r1_40um_d0.05.tiff'));
I = imcrop(I, [2042-hw, 2018-hw, 2*hw, 2*hw]);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, 'A_dwell_8us.tiff'));

I = imread(fullfile(working_dir, 'def_s16_r1_40um_d0.05.tiff'));
I = imcrop(I, [2052-hw, 2010-hw, 2*hw, 2*hw]);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, 'A_dwell_16us.tiff'));

%% [B] Reducing FOV / increasing mag can lead to big improvement
close all;

I = imread(fullfile(working_dir, 'ref_s4_r1_40um_d0.16.tiff'));
I = imcrop(I, [1885-hw, 2050-hw, 2*hw, 2*hw]);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, 'B_fov_40um.tiff'));

I = imread(fullfile(working_dir, 'def_s4_r1_20um_d0.05.tiff'));
I = imcrop(I, [1040-hw, 977-hw, 2*hw, 2*hw]);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, 'B_fov_20um.tiff'));

I = imread(fullfile(working_dir, 'def_s4_r1_10um_d0.05.tiff'));
I = imcrop(I, [531-hw, 463-hw, 2*hw, 2*hw]);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, 'B_fov_10um.tiff'));

%% [C] Snake vs Raster are similar, if total time is similar
close all;

I = imread(fullfile(working_dir, 'def_s8_r1_40um_d0.05.tiff'));
I = imcrop(I, [2042-hw, 2018-hw, 2*hw, 2*hw]);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, 'C_raster_8us.tiff'));

I = imread(fullfile(working_dir, 'ref_s4_r0_40um_d0.16.tiff'));
I = imcrop(I, [1893-hw, 2041-hw, 2*hw, 2*hw]);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, 'C_snake_4x2us.tiff'));

%% [0] Built-in for reference
close all;

I = imread(fullfile(working_dir, 'ref_s4_r1_40um_d0.16.tiff'));
I = imcrop(I, [1885-hw, 2050-hw, 2*hw, 2*hw]);
figure;
imshow(I);

I = imread(fullfile(working_dir, 'built-in image 0830.tif'));
I = imcrop(I, [1855-hw, 2070-hw, 2*hw, 2*hw]);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, '0_builtin_ref.tiff'));
%%
close all;


