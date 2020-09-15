addChenFunction;
%% (1) Imaged by Python function, dwell time = 1us, total time about 16 sec
close all;

cd('D:\p\python\PWD\Tescan-SharkSEM\img_20200914\compare image quality\auto_dwell_1us')

d1 = load('cali_r0c0_A.mat')
d2 = load('cali_r0c1.mat')

u1 = nanmean(d1.u(:));
std1 = nanstd(d1.exx(:));
myplot(d1.x, d1.y, d1.exx);
caxis([-0.02 0.02]);
title_str = ['exx, u = ',num2str(u1,'%3.0f'),' pxls, std = ',num2str(std1,'%3.4f')];
title(title_str, 'fontweight', 'normal');

u2 = nanmean(d2.u(:));
std2 = nanstd(d2.exx(:));
myplot(d2.x, d2.y, d2.exx);
caxis([-0.02 0.02]);
title_str = ['exx, u = ',num2str(u2,'%3.0f'),' pxls, std = ',num2str(std2,'%3.4f')];
title(title_str, 'fontweight', 'normal');

%% (2) Imaged by built-in function, scan speed = 3, total time about 18 sec
cd('D:\p\python\PWD\Tescan-SharkSEM\img_20200914\compare image quality\speed_3')

d1 = load('image_0002.mat')
d2 = load('image_0005.mat')

u1 = nanmean(d1.u(:));
std1 = nanstd(d1.exx(:));
myplot(d1.x, d1.y, d1.exx);
caxis([-0.02 0.02]);
title_str = ['exx, u = ',num2str(u1,'%3.0f'),' pxls, std = ',num2str(std1,'%3.4f')];
title(title_str, 'fontweight', 'normal');

u2 = nanmean(d2.u(:));
std2 = nanstd(d2.exx(:));
myplot(d2.x, d2.y, d2.exx);
caxis([-0.02 0.02]);
title_str = ['exx, u = ',num2str(u2,'%3.0f'),' pxls, std = ',num2str(std2,'%3.4f')];
title(title_str, 'fontweight', 'normal');

%% (3) Imaged by built-in function, scan speed = 3, line integration = 5, total time about 90 sec
cd('D:\p\python\PWD\Tescan-SharkSEM\img_20200914\compare image quality\speed_3 line_int_5')

d1 = load('image_0004.mat')
d2 = load('image_0007.mat')

u1 = nanmean(d1.u(:));
std1 = nanstd(d1.exx(:));
myplot(d1.x, d1.y, d1.exx);
caxis([-0.02 0.02]);
title_str = ['exx, u = ',num2str(u1,'%3.0f'),' pxls, std = ',num2str(std1,'%3.4f')];
title(title_str, 'fontweight', 'normal');

u2 = nanmean(d2.u(:));
std2 = nanstd(d2.exx(:));
myplot(d2.x, d2.y, d2.exx);
caxis([-0.02 0.02]);
title_str = ['exx, u = ',num2str(u2,'%3.0f'),' pxls, std = ',num2str(std2,'%3.4f')];
title(title_str, 'fontweight', 'normal');



