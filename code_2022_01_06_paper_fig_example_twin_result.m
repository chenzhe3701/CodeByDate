% [Fig 4]
% example of EBSD parent grain, EBSD child grain, grain-level twin map,
% pixel-level twin map.
% sample Mg4Al_C3, iE=2

addChenFunction;
clear; clc; close all;
output_dir = 'C:\Users\chenz\Work\Data\0_temp_output\fig 4';
mkdir(output_dir);
%%

data_dir = 'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis';
variant_dir = 'H:\Other computers\My Laptop w541 2022\zhec umich Drive\All twin variant maps cleaned'; 

id_p = 101;
id_c = [294,395,398,399,400];

indR_min = 280;
indR_max = 400;
indC_min = 260;
indC_max = 380;

d = matfile(fullfile(variant_dir, 'Mg4Al_C3_variant_maps.mat'));
vp = d.variant_point_wise;
vg = d.variant_grain_wise;
vp = vp{3};
vg = vg{3};
vp = vp(indR_min:indR_max, indC_min:indC_max);
vg = vg(indR_min:indR_max, indC_min:indC_max);

d = matfile(fullfile(data_dir, 'Mg4Al_C3_parent_grain_file_iE_2.mat'));
ID = d.ID;
phi1 = d.phi1;
phi = d.phi;
phi2 = d.phi2;
ID = ID(indR_min:indR_max, indC_min:indC_max);
phi1 = phi1(indR_min:indR_max, indC_min:indC_max);
phi = phi(indR_min:indR_max, indC_min:indC_max);
phi2 = phi2(indR_min:indR_max, indC_min:indC_max);
boundary = find_one_boundary_from_ID_matrix(ID);

d = matfile(fullfile(data_dir, 'Mg4Al_C3_grain_file_iE_2.mat')); 
ID_c = d.ID;
ID_c = ID_c(indR_min:indR_max, indC_min:indC_max);
boundary_c = find_one_boundary_from_ID_matrix(ID_c);

% make IPF
[nR,nC] = size(ID);
for iR = 1:nR
   for iC = 1:nC
        RGB(iR,iC,:) = reshape(calculate_IPF_color_hcp([phi1(iR,iC), phi(iR,iC), phi2(iR,iC)], [0 0 0], [0 0 1]), 1, 1, 3);
   end
end

IPF_p = RGB;
for iR = 1:nR
   for iC = 1:nC
       if boundary(iR,iC) == 1
        IPF_p(iR,iC,:) = reshape([1 1 1], 1, 1, 3);
       end
   end
end

IPF_c = RGB;
for iR = 1:nR
   for iC = 1:nC
       if boundary_c(iR,iC) == 1
        IPF_c(iR,iC,:) = reshape([1 1 1], 1, 1, 3);
       end
   end
end

% figure;
% imshow(IPF_p);
% figure;
% imshow(IPF_c)
% 
% myplot(vg, boundary);
% myplot(vp, boundary);

%%
close all;

% IPF parent
figure; imagesc(IPF_p);
axis equal; axis off
[x,y] = meshgrid((1:size(ID,2))-1, (1:size(ID,1))-1);
label_map_with_ID(x,y,ID,gcf,101,'k',16,1)

print(fullfile(output_dir, 'fig 4 parent IPF.tif'),'-dtiff');


% IPF child
figure; imagesc(IPF_c);
axis equal; axis off
[x,y] = meshgrid((1:size(ID,2))-1, (1:size(ID,1))-1);
label_map_with_ID(x,y,ID_c,gcf,[394,395,398,399,400],'k',16,1)

print(fullfile(output_dir, 'fig 4 child IPF.tif'),'-dtiff');


% variant grain level
close all;
map = vg;
map(boundary==1) = 7;

figure; 
imagesc(map); 
axis equal; axis off;
set(gca,'fontsize',18);

cmap = parula(7);
cmap(1,:) = [0.25 0.25 0.25];
cmap(8,:) = [1 1 1];
colormap(cmap);

cbar = colorbar;
caxis([0 7+1]);
set(cbar,'limit',[1 7],'Ticks',[1.5:6.5],'TickLabels',{'1','2','3','4','5','6'});

print(fullfile(output_dir, 'fig 4 variant grain.tif'),'-dtiff');


% variant pixel level
close all;
map = vp;
map(boundary==1) = 7;

figure; 
imagesc(map); 
axis equal; axis off;
set(gca,'fontsize',18);

cmap = parula(7);
cmap(1,:) = [0.25 0.25 0.25];
cmap(8,:) = [1 1 1];
colormap(cmap);

cbar = colorbar;
caxis([0 7+1]);
set(cbar,'limit',[1 7],'Ticks',[1.5:6.5],'TickLabels',{'1','2','3','4','5','6'});

print(fullfile(output_dir, 'fig 4 variant pixel.tif'),'-dtiff');

close all;

%% crop
if 1
    I = imread(fullfile(output_dir, 'fig 4 parent IPF.tif'));
    I = imcrop(I, [185 49 538 538]);
    imwrite(I, fullfile(output_dir, 'fig 4 parent IPF.tif'));
    
    I = imread(fullfile(output_dir, 'fig 4 child IPF.tif'));
    I = imcrop(I, [185 49 538 538]);
    imwrite(I, fullfile(output_dir, 'fig 4 child IPF.tif'));
    
    I = imread(fullfile(output_dir, 'fig 4 variant grain.tif'));
    I = imcrop(I, [165 76 633 482]);
    imwrite(I, fullfile(output_dir, 'fig 4 variant grain.tif'));
    
    I = imread(fullfile(output_dir, 'fig 4 variant pixel.tif'));
    I = imcrop(I, [165 76 633 482]);
    imwrite(I, fullfile(output_dir, 'fig 4 variant pixel.tif')); 
end