% output EBSD based data, and try DIC

output_dir = 'C:\Users\ZheChen\Desktop\try dic';
mkdir(output_dir);

%% undeformed EBSD
iE_0 = 0;
d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE_0),'.mat']));
ID = d.ID;    % parent grain
X = d.x;
Y = d.y;
boundary_0 = find_one_boundary_from_ID_matrix(ID);

m = city_block(boundary_0);
I = mat_to_image(m, [0 70], 'uint8');

myplot(boundary_0);
figure; imagesc(I);

fname = ['fig iE_',num2str(iE_0)];
imwrite(I, fullfile(output_dir, [fname,'.tiff']),'tiff');
close all;

%% deformed EBSD
iE = 3;
d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
ID = d.ID;    % parent grain
X = d.x;
Y = d.y;
boundary = find_one_boundary_from_ID_matrix(ID);

m = city_block(boundary);
I = mat_to_image(m, [0 70], 'uint8');

myplot(boundary);
figure; imagesc(I);

fname = ['fig iE_',num2str(iE)];
imwrite(I, fullfile(output_dir, [fname,'.tiff']),'tiff');
close all;

%% after doing DIC
close all;

iE = 3;
fname = ['fig iE_',num2str(iE)];
d = matfile(fullfile(output_dir, [fname,'.mat']));
x = d.x;
y = d.y;
u = d.u;
v = d.v;
sigma = d.sigma;
% u(sigma==-1) = nan;
% v(sigma==-1) = nan;

% data points on the undeformed map (boundary_0) have coordinates from the deformed map (xp, yp)  
xp = x+u;   
yp = y+v;

[~, locb_x] = ismember(x(1,:), X(1,:));
[~, locb_y] = ismember(y(:,1), Y(:,1));
boundary_0_cropped = boundary_0(locb_y, locb_x);    % crop

% interp data to grid coordinate
% boundary_0_interpped = interp_data(xp, yp, boundary_0_cropped, X,Y, [], 'fit', 'nearest');
F = scatteredInterpolant(xp(:), yp(:), boundary_0_cropped(:));
boundary_0_interpped = F(X, Y);

figure; hold on;
% (1) EBSD measurement at iE_0
% ind = (boundary_0==1);
% plot(X(ind), Y(ind), '.k', 'markersize', 3);  

% (2) EBSD measurement at iE
ind = (boundary==1);
plot(X(ind), Y(ind), '.r', 'markersize', 3);  

% (3) plot boundary at iE_0 using DIC determined coordinate at iE
ind = (boundary_0_cropped==1);
plot(xp(ind), yp(ind), '.b', 'markersize', 3);

% (4) interp DIC measurement to grid coordinate
% ind = (boundary_0_interpped==1);
% plot(X(ind), Y(ind), '.k', 'markersize', 30);

% legend('EBSD measurement (ref)', ...
%     'EBSD measurement (deformed)', ...
%     'ref data with deformed position by DIC', ...
%     'ref data with deformed position by DIC, to grids', ...
%     'location', 'northoutside'); 
axis equal;








