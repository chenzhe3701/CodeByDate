% Make fig 8 using UM134_Mg_C3

fig_dir = 'C:\Users\chenz\OneDrive\Documents\Grain Size and Alloying Effect Paper\fig 8 combined maps\um134 c3';

for iE = 0:13
    f1 = ['variant_pt_wise_iE=',num2str(iE),'.tif'];
    f2 = ['UM134_Mg_C3 frd twinned iE=',num2str(iE),'.tiff'];
    f3 = ['UM134_Mg_C3 grain label iE=',num2str(iE),'.tiff'];
    
    I1 = imread(fullfile(fig_dir,f1));
    I2 = imread(fullfile(fig_dir,f2));
    I3 = imread(fullfile(fig_dir,f3));



    I1 = imcrop(I1, [201,71,765,765]);
    I2 = imcrop(I2, [379,71,765,765]);
    I3 = imcrop(I3, [541,71,765,765]);

    I1([1,end],:,:) = 0;
    I1(:,[1,end],:) = 0;
    I2([1,end],:,:) = 0;
    I2(:,[1,end],:) = 0;
    I3([1,end],:,:) = 0;
    I3(:,[1,end],:) = 0;

    I = uint8(255*ones(785,765*3+40,3));
    I(10:775, 10:775, :) = I1;
    I(10:775, 785:1550, :) = I2;
    I(10:775, 1560:2325, :) = I3;

    imwrite(I, fullfile(fig_dir, ['img_',num2str(iE),'.tiff']));
end
