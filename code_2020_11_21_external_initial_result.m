
%% (1) For the external scan trial test on 2020-11-18, the images with external scan show streaks in the horizontal direction.

working_dir = 'E:\zhec umich Drive\2020-11-18 Ext vs Builtin on Tescan\Compare AOI_1\';
cd(working_dir);

image_name = 'ext_image_001 raw.tiff';
I_r = imread(fullfile(working_dir, image_name));

image_name = 'ext_image_001.tiff';
I = imread(fullfile(working_dir, image_name));

ref_name = 'built-in image-001.tif';
R = imread(fullfile(working_dir, ref_name));

%% image I_r is auto B&C adjusted by photoshop. compare their value
figure;
histogram(I_r(:));
figure;
histogram(I(:));

%% try initial filter
close all;
mu = [0, 0];
sigma = [0.5,0; 0,5]*1;

x = -1:1;
y = -5:5;
[X, Y] = meshgrid(x, y);
XX = [X(:), Y(:)];
YY = mvnpdf(XX, mu, sigma); % multivariat normal distribution

initial_psf = reshape(YY, length(y), length(x))
II = imfilter(R, initial_psf);

% Image by built-in scanner
figure;
imshow(R);
set(gca,'xlim',[2320, 2620], 'ylim', [3100, 3400]);
title('built-in scanner');

% Image by external scanner
figure;
imshow(I);
set(gca,'xlim',[2385, 2685], 'ylim', [3155, 3455]);
title('ext scanner');

% forward, image by built-in scanner processed by a 2-d filter corresponding to multivariant normal distribution
figure;
imshow(auto_BC(II));
set(gca,'xlim',[2320, 2620], 'ylim', [3100, 3400]);
title('built-in + filter');

% The filter
figure;imagesc(initial_psf);

%% We can try to solve the filter, and the deconvoluted image. Currently the result is not good.
[J P] = deconvblind(I, initial_psf, 2);

figure;
imshow(J);
title('deconvoluted');

%% (2.1) Look at DIC results for sample translation results, built-in scanner
built_in_scanner_dir = 'E:\zhec umich Drive\2020-11-18 Ext vs Builtin on Tescan\Ref crosshatch\';
for ii = 2:5
    close all;
    built_in_image_name = ['built-in image-00',num2str(ii),'.mat'];
    
    data = matfile(fullfile(built_in_scanner_dir,built_in_image_name));
    u = data.u;
    v = data.v;
    exx = data.exx;
    exy = data.exy;
    eyy = data.eyy;
    sigma = data.sigma;
    
    ind = sigma == -1;
    u(ind) = nan;
    v(ind) = nan;
    exx(ind) = nan;
    exy(ind) = nan;
    eyy(ind) = nan;
    
    myplot(u);
    set(gca,'fontsize',18,'xTick',[],'yTick',[]);
    print(fullfile(built_in_scanner_dir,['u_00',num2str(ii),'.tiff']),'-dtiff');
    
    myplot(v);
    set(gca,'fontsize',18,'xTick',[],'yTick',[]);
    print(fullfile(built_in_scanner_dir,['v_00',num2str(ii),'.tiff']),'-dtiff');
    
    myplot(exx);
    caxis([-0.02, 0.02]);
    set(gca,'fontsize',18,'xTick',[],'yTick',[]);
    print(fullfile(built_in_scanner_dir,['exx_00',num2str(ii),'.tiff']),'-dtiff');
    
    myplot(exy);
    caxis([-0.02, 0.02]);
    set(gca,'fontsize',18,'xTick',[],'yTick',[]);
    print(fullfile(built_in_scanner_dir,['exy_00',num2str(ii),'.tiff']),'-dtiff');
    
    myplot(eyy);
    caxis([-0.02, 0.02]);
    set(gca,'fontsize',18,'xTick',[],'yTick',[]);
    print(fullfile(built_in_scanner_dir,['eyy_00',num2str(ii),'.tiff']),'-dtiff');
end
close all;

% (2.2) Look at DIC results for sample translation results, external scanner
external_scanner_dir = 'E:\zhec umich Drive\2020-11-18 Ext vs Builtin on Tescan\Ext scan image\';
for ii = 2:5
    close all;
    external_scanner_image_name = ['ext_image_00',num2str(ii),'.mat'];
    
    data = matfile(fullfile(external_scanner_dir, external_scanner_image_name));
    u = data.u;
    v = data.v;
    exx = data.exx;
    exy = data.exy;
    eyy = data.eyy;
    sigma = data.sigma;
    
    ind = sigma == -1;
    u(ind) = nan;
    v(ind) = nan;
    exx(ind) = nan;
    exy(ind) = nan;
    eyy(ind) = nan;
    
    myplot(u);
    myplot(v);
    myplot(exx);
    myplot(exy);
    myplot(eyy);
    
    myplot(u);
    set(gca,'fontsize',18,'xTick',[],'yTick',[]);
    print(fullfile(external_scanner_dir,['u_00',num2str(ii),'.tiff']),'-dtiff');
    
    myplot(v);
    set(gca,'fontsize',18,'xTick',[],'yTick',[]);
    print(fullfile(external_scanner_dir,['v_00',num2str(ii),'.tiff']),'-dtiff');
    
    myplot(exx);
    set(gca,'fontsize',18,'xTick',[],'yTick',[]);
    caxis([-0.02, 0.02]);
    print(fullfile(external_scanner_dir,['exx_00',num2str(ii),'.tiff']),'-dtiff');
    
    myplot(exy);
    set(gca,'fontsize',18,'xTick',[],'yTick',[]);
    caxis([-0.02, 0.02]);
    print(fullfile(external_scanner_dir,['exy_00',num2str(ii),'.tiff']),'-dtiff');
    
    myplot(eyy);
    set(gca,'fontsize',18,'xTick',[],'yTick',[]);
    caxis([-0.02, 0.02]);
    print(fullfile(external_scanner_dir,['eyy_00',num2str(ii),'.tiff']),'-dtiff');
end
close all;










