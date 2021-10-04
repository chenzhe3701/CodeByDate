% find control points for Mohsen's Mg4Al
close all;
clc;

% EBSD data
ebsd_dir = 'E:\zhec umich Drive\2021-10-01 Mg4Al_m1 insitu SEM-DIC\EBSD Data';

f1 = fullfile(ebsd_dir, 'Mg4Al_grain_file_type_1.txt');
f2 = fullfile(ebsd_dir, 'Mg4Al_grain_file_type_2.txt');

d = grain_file_to_data(f1,f2);

ID = d.ID;
boundary = find_one_boundary_from_ID_matrix(ID);

figure;
imagesc(boundary);
[f1,a1,c1] = myplot(boundary);

%% DIC data
dic_data_dir = 'E:\zhec umich Drive\2021-10-01 Mg4Al_m1 insitu SEM-DIC\SEM Data\stitched DIC';
d = matfile(fullfile(dic_data_dir, 'DIC_merged_3.mat'));
x = d.x;
y = d.y;
exx = d.exx;
exy = d.exy;
eyy = d.eyy;
sigma = d.sigma;

exx(sigma==-1) = nan;
exy(sigma==-1) = nan;
eyy(sigma==-1) = nan;

eeff = calculate_effective_strain(exx,exy,eyy);

[f2,a2,c2] = myplot(x,y,eeff);
% myplot(x,y,exx);
% myplot(x,y,exy);
% myplot(x,y,eyy);
%% SEM image
% sem_dir = 'E:\zhec umich Drive\2021-10-01 MgAl insitu SEM-DIC\SEM Data\stitched images';
% I = imread(fullfile(sem_dir, 'img_full_res_0.tif'));
% figure;
% imshow(I);

%% record control points
cpSEM = [1812, 2652;
    15370, 2322;
    5872, 9902;
    16600, 11370];

cpEBSD = [28, 82;
    231, 81;
    89, 188;
    247, 217];

plot(a1, cpEBSD(:,1), cpEBSD(:,2), 'or');
plot(a2, cpSEM(:,1), cpSEM(:,2), 'or');
