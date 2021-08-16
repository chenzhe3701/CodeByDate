% try to crop a relatively small image for Mg4Al_S1 sample to try the ALDIC
% codes

working_dir = 'C:\Users\ZheChen\Desktop\try dic';
cd(working_dir);

%%
I = imread('Mg4Al_S1_s0_r0c1.tif');
Ic = I(2101:3200, 501:1600);
imwrite(Ic, 'Mg4Al_img_a.tif');

I = imread('Mg4Al_S1_s2_r0c1.tif');
Ic = I(2101:3200, 501:1600);
imwrite(Ic, 'Mg4Al_img_b.tif');