
working_dir = 'E:\zhec umich Drive\2021-05-25 test image quality';
dic_file = '15kv_11bi_15wd_speed3_line5_fov45.mat';
addChenFunction;

d = matfile(fullfile(working_dir, dic_file));
x = d.x;
y = d.y;
exx = d.exx;
u = d.u;
sigma = d.sigma;

exx(sigma==-1) = nan;
u(sigma==-1) = nan;

myplot(x,y,exx);
print(fullfile(working_dir, 'exx.tiff'),'-dtiff');

myplot(x,y,u);
print(fullfile(working_dir, 'u.tiff'),'-dtiff');

close all;