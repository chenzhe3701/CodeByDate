%% 2021-08-29
% analyze the quality of external scan after moving to Dow
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-08-27 try external in Dow';
output_dir = fullfile(working_dir, 'output');
mkdir(output_dir);

%% (1) plot the distortion of the builtin image 103
close all;
d = matfile(fullfile(working_dir, 'built-in image_103.mat'));
x = d.x;
y = d.y;
u = d.u;
v = d.v;
exx = d.exx;
sigma = d.sigma;

x(sigma==-1) = nan;
y(sigma==-1) = nan;
u(sigma==-1) = nan;
v(sigma==-1) = nan;
exx(sigma==-1) = nan;

myplot(x,y,u);
set(gca,'fontsize',16);xlabel('x (pixel)');ylabel('y (pixel)');
print(fullfile(output_dir, 'built-in image 103 u.tiff'),'-dtiff');

myplot(x,y,v);
set(gca,'fontsize',16);xlabel('x (pixel)');ylabel('y (pixel)');
print(fullfile(output_dir, 'built-in image 103 v.tiff'),'-dtiff');

myplot(x,y,exx);
set(gca,'fontsize',16);xlabel('x (pixel)');ylabel('y (pixel)');
print(fullfile(output_dir, 'built-in image 103 exx.tiff'),'-dtiff');

%% (1 additional) exagerate displacement, illustrate the deformed position of a regular grid
close all;
indR_min = 1;
indR_max = 70;
indC_min = 1;
indC_max = 70;
x_local = x(indR_min:indR_max, indC_min:indC_max);
y_local = y(indR_min:indR_max, indC_min:indC_max);
u_local = u(indR_min:indR_max, indC_min:indC_max);
v_local = v(indR_min:indR_max, indC_min:indC_max);

xp = x_local + 15 * u_local;
yp = y_local + 15 * v_local;

[nR,nC] = size(x_local);

figure;
hold on;
for iR = 1:1:nR
    plot(xp(iR,:), yp(iR,:), '.-k');
end
for iC = 1:1:nC
   plot(xp(:,iC), yp(:,iC), '.-k'); 
end
set(gca,'ydir','reverse');
axis equal;

print(fullfile(output_dir, 'exaggerate displacement 15x builtin image 103.tiff'), '-dtiff');

%% (2) plot the distortion of the external scan image 103
close all;
d = matfile(fullfile(working_dir, 'external_image_103.mat'));
x = d.x;
y = d.y;
u = d.u;
v = d.v;
exx = d.exx;
sigma = d.sigma;

x(sigma==-1) = nan;
y(sigma==-1) = nan;
u(sigma==-1) = nan;
v(sigma==-1) = nan;
exx(sigma==-1) = nan;

myplot(x,y,u);
set(gca,'fontsize',16);xlabel('x (pixel)');ylabel('y (pixel)');
print(fullfile(output_dir, 'external image 103 u.tiff'),'-dtiff');

myplot(x,y,v);
set(gca,'fontsize',16);xlabel('x (pixel)');ylabel('y (pixel)');
print(fullfile(output_dir, 'external image 103 v.tiff'),'-dtiff');

myplot(x,y,exx);
set(gca,'fontsize',16);xlabel('x (pixel)');ylabel('y (pixel)');
print(fullfile(output_dir, 'external image 103 exx.tiff'),'-dtiff');

%% (3) show the actual image of a small area, to compare
close all;
% (3.1) external scan image
I = imread(fullfile(working_dir,'external_image_101.tiff'));
x_target = 2048;
y_target = 2048;
hw = 100;
I = imcrop(I, [x_target-hw, y_target-hw, 2*hw, 2*hw]);
figure; imshow(I);
imwrite(I, fullfile(output_dir, 'f1 ext fov=40um dwell=4.tiff'));

% (3.2) stock scan image
d = matfile(fullfile(working_dir, 'built-in image_101.mat'));
x = d.x;
y = d.y;
u = d.u;
v = d.v;
% find where the feature located on the stock scan image, by looking its
% u/v value wrt the external scan image
[~,indc] = min(abs(x(1,:) - x_target)); 
[~,indr] = min(abs(y(:,1) - y_target));
x_tgt = x_target + u(indr,indc);
y_tgt = y_target + v(indr,indc);
I = imread(fullfile(working_dir, 'built-in image_101.tif'));
I = imcrop(I, [x_tgt-hw, y_tgt-hw, 2*hw, 2*hw]);
figure; imshow(I);
imwrite(I, fullfile(output_dir, 'f0 builtin fov=40um.tiff'));

% (3.3) external image, long dwell time = 10
d = matfile(fullfile(working_dir, 'external_image_s10_101.mat'));
x = d.x;
y = d.y;
u = d.u;
v = d.v;
% find where the feature located on the stock scan image, by looking its
% u/v value wrt the external scan image
[~,indc] = min(abs(x(1,:) - x_target)); 
[~,indr] = min(abs(y(:,1) - y_target));
x_tgt = x_target + u(indr,indc);
y_tgt = y_target + v(indr,indc);
I = imread(fullfile(working_dir, 'external_image_s10_101.tiff'));
I = imcrop(I, [x_tgt-hw, y_tgt-hw, 2*hw, 2*hw]);
figure; imshow(I);
imwrite(I, fullfile(output_dir, 'f2 ext fov=40um dwell=10.tiff'));

% (3.4) external image, long dwell time = 20
d = matfile(fullfile(working_dir, 'external_image_s20_101.mat'));
x = d.x;
y = d.y;
u = d.u;
v = d.v;
% find where the feature located on the stock scan image, by looking its
% u/v value wrt the external scan image
[~,indc] = min(abs(x(1,:) - x_target)); 
[~,indr] = min(abs(y(:,1) - y_target));
x_tgt = x_target + u(indr,indc);
y_tgt = y_target + v(indr,indc);
I = imread(fullfile(working_dir, 'external_image_s20_101.tiff'));
I = imcrop(I, [x_tgt-hw, y_tgt-hw, 2*hw, 2*hw]);
figure; imshow(I);
imwrite(I, fullfile(output_dir, 'f3 ext fov=40um dwell=20.tiff'));

% (3.5) external image, regular dwell = 4, but small area of 10um
I = imread(fullfile(working_dir,'external_10um_101.tiff'));
x_tgt = 547;
y_tgt = 514;
hw = 100;
I = imcrop(I, [x_tgt-hw, y_tgt-hw, 2*hw, 2*hw]);
figure; imshow(I);
imwrite(I, fullfile(output_dir, 'f4 ext fov=10um dwell=4.tiff'));

%%
close all;

