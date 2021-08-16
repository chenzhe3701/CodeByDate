%% data for DIC with boundary. Collaboration with Jin Yang
clear; clc;
addChenFunction;
data_dir = 'E:\Ti7Al_B6_insitu_tension';
output_dir = 'E:\zhec umich Drive\0_temp_output\Ti7Al_B6_develop_dic';
mkdir(output_dir);

%% reconstructed EBSD data, has the same resolution and coordinate as stitched DIC data
ebsd_data = matfile(fullfile(data_dir, 'Analysis_by_Matlab/Ti7Al_B6_EbsdToSemForTraceAnalysis_GbAdjusted.mat'));
X = ebsd_data.X;
Y = ebsd_data.Y;
boundaryTF = ebsd_data.boundaryTF;
boundaryTFB = ebsd_data.boundaryTFB;
ID = ebsd_data.ID;

%% transX and transY for stitching different FOV of SEM images
transXY_data = matfile(fullfile(data_dir, 'SEM Data', 'stitch gage image', 'translations_searched_vertical_stop_0.mat'));
transX = transXY_data.transX;
transY = transXY_data.transY;
% this is how transX/Y pre-processed in my code
dic_step = 2;
transX = dic_step * round(transX/dic_step);
transY = dic_step * round(transY/dic_step);

%% illustrate DIC data of selected iE, and the 36 FOVs
iE = 11;
uvTrans_data = fullfile(data_dir, 'SEM Data\stitched_DIC\uvTrans', ['uvTrans_e',num2str(iE),'.mat']);
load(uvTrans_data, 'uTrans', 'vTrans'); % this is how the 'u' and 'v' should be corrected

% stitched dic data
dic_data = matfile(fullfile(data_dir, 'SEM Data\stitched_DIC', ['_',num2str(iE),'.mat']));
exx = dic_data.exx;

% illustrate all FOVs
myplot(X,Y,exx,boundaryTFB);
for iR = 1:6
    for iC = 1:6
        drawrectangle('Position', [transX(iR,iC), transY(iR,iC), 4096, 4096], 'color', 'r','InteractionsAllowed','none');
    end
end

print(fullfile(output_dir,['exx_all_fovs_iE=',num2str(iE),'.tiff']),'-dtiff');
close;
%% [fov 1:] Select a FOV[iR,iC] = [1,1] 
iR = 1;
iC = 1;
B = 1;

mkdir(fullfile(output_dir,['r',num2str(iR),'c',num2str(iC)]));

% for each FOV image, coordinates are 0:4095, 0:4095
[x_SEM_FOV, y_SEM_FOV] = meshgrid(0:4095, 0:4095);

xyi = 1;    % for this sample, each FOV has coordinates [7:2:4087, 7:2:4087] 

% crop the boundary/EBSD data within this area for the selected FOV[iR,iC]
X_crop_min = 0 + transX(iR+B,iC+B);  
Y_crop_min = 0 + transY(iR+B,iC+B);
X_crop_max = X_crop_min + 4095;
Y_crop_max = Y_crop_min + 4095;

% indices to crop
indC_min = find( X(1,:)>=X_crop_min, 1, 'first');
indC_max = find( X(1,:)<=X_crop_max, 1, 'last');
indR_min = find( Y(:,1)>=Y_crop_min, 1, 'first');
indR_max = find( Y(:,1)<=Y_crop_max, 1, 'last');
            
% crop   
X_crop = X(indR_min:indR_max, indC_min:indC_max);
X_crop = X_crop - X_crop(1) + xyi;
Y_crop = Y(indR_min:indR_max, indC_min:indC_max);
Y_crop = Y_crop - Y_crop(1) + xyi;

boundary_crop = boundaryTF(indR_min:indR_max, indC_min:indC_max);
ID_crop = ID(indR_min:indR_max, indC_min:indC_max);

% boundary with the same resolution as the reference image
boundary = interp_data(X_crop, Y_crop, boundary_crop, x_SEM_FOV, y_SEM_FOV, [], 'interp', 'nearest');
boundary(isnan(boundary)) = 0;

save(fullfile(output_dir, ['Ti7Al_B6_boundary_r',num2str(iR),'c',num2str(iC),'.mat']), 'boundary');

myplot(x_SEM_FOV, y_SEM_FOV, boundary);
print(fullfile(output_dir, ['boundary_r',num2str(iR),'c',num2str(iC),'.tiff']),'-dtiff');
close;

% output for each iE, 2041x2041: x,y,u,v,exx,exy,eyy,sigma
for iE = 0:14    
    fov_data = load(fullfile(data_dir, 'SEM Data\byFov_IMAGE_DELETED_TO_REDUCE_SIZE', ['r',num2str(iR),'c',num2str(iC)], ...
        ['Ti7Al_B6_e',num2str(iE),'_r',num2str(iR),'c',num2str(iC),'.mat']));
    x = fov_data.x;
    y = fov_data.y;
    u = fov_data.u;
    v = fov_data.v;
    exx = fov_data.exx;
    exy = fov_data.exy;
    eyy = fov_data.eyy;
    sigma = fov_data.sigma;
    save(fullfile(output_dir, ['r',num2str(iR),'c',num2str(iC)], ['Ti7Al_B6_e',num2str(iE),'_r',num2str(iR),'c',num2str(iC),'.mat']), ...
        'x','y','u','v','exx','exy','eyy','sigma');
end

%% [fov 2:] Select a FOV[iR,iC] = [1,1], [1,2]
iR = 1;
iC = 2;
B = 1;

mkdir(fullfile(output_dir,['r',num2str(iR),'c',num2str(iC)]));

% for each FOV image, coordinates are 0:4095, 0:4095
[x_SEM_FOV, y_SEM_FOV] = meshgrid(0:4095, 0:4095);

xyi = 1;    % for this sample, each FOV has coordinates [7:2:4087, 7:2:4087] 

% crop the boundary/EBSD data within this area for the selected FOV[iR,iC]
X_crop_min = 0 + transX(iR+B,iC+B);  
Y_crop_min = 0 + transY(iR+B,iC+B);
X_crop_max = X_crop_min + 4095;
Y_crop_max = Y_crop_min + 4095;

% indices to crop
indC_min = find( X(1,:)>=X_crop_min, 1, 'first');
indC_max = find( X(1,:)<=X_crop_max, 1, 'last');
indR_min = find( Y(:,1)>=Y_crop_min, 1, 'first');
indR_max = find( Y(:,1)<=Y_crop_max, 1, 'last');
            
% crop   
X_crop = X(indR_min:indR_max, indC_min:indC_max);
X_crop = X_crop - X_crop(1) + xyi;
Y_crop = Y(indR_min:indR_max, indC_min:indC_max);
Y_crop = Y_crop - Y_crop(1) + xyi;

boundary_crop = boundaryTF(indR_min:indR_max, indC_min:indC_max);
ID_crop = ID(indR_min:indR_max, indC_min:indC_max);

% boundary with the same resolution as the reference image
boundary = interp_data(X_crop, Y_crop, boundary_crop, x_SEM_FOV, y_SEM_FOV, [], 'interp', 'nearest');
boundary(isnan(boundary)) = 0;

save(fullfile(output_dir, ['Ti7Al_B6_boundary_r',num2str(iR),'c',num2str(iC),'.mat']), 'boundary');

myplot(x_SEM_FOV, y_SEM_FOV, boundary);
print(fullfile(output_dir, ['boundary_r',num2str(iR),'c',num2str(iC),'.tiff']),'-dtiff');
close;

% output for each iE, 2041x2041: x,y,u,v,exx,exy,eyy,sigma
for iE = 0:14    
    fov_data = load(fullfile(data_dir, 'SEM Data\byFov_IMAGE_DELETED_TO_REDUCE_SIZE', ['r',num2str(iR),'c',num2str(iC)], ...
        ['Ti7Al_B6_e',num2str(iE),'_r',num2str(iR),'c',num2str(iC),'.mat']));
    x = fov_data.x;
    y = fov_data.y;
    u = fov_data.u;
    v = fov_data.v;
    exx = fov_data.exx;
    exy = fov_data.exy;
    eyy = fov_data.eyy;
    sigma = fov_data.sigma;
    save(fullfile(output_dir, ['r',num2str(iR),'c',num2str(iC)], ['Ti7Al_B6_e',num2str(iE),'_r',num2str(iR),'c',num2str(iC),'.mat']), ...
        'x','y','u','v','exx','exy','eyy','sigma');
end

%% OK to copy images manually

