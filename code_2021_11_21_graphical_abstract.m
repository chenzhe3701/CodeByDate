% make strain + twin map for graphical abstract

clear;
addChenFunction;
output_dir = 'E:\zhec umich Drive\0_temp_output';

% looks like have to include this part to read the sample name.
load_settings('D:\p\m\DIC_Analysis\setting_for_real_samples\WE43_T6_C1_setting.mat', ...
    'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor','strainPauses');

% load previous data and settings
load('D:\WE43_T6_C1\Analysis_by_Matlab_after_realign\WE43_T6_C1_EbsdToSemForTraceAnalysis_GbAdjusted', ...
    'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');

load('D:\WE43_T6_C1\Analysis_2021_09\WE43_T6_C1_new_variant_map.mat', ...
    'struCell','variantMapCleanedCell','trueTwinMapCell');

%%
variant = variantMapCleanedCell{4};
exx = matfile('D:\WE43_T6_C1\SEM Data\stitched_DIC\_4_v73.mat').exx;

indr = 3200:3900;
indc = 9300:10000;

variant_local = variant(indr,indc);
exx_local = exx(indr,indc);
boundary_local = boundaryTFB(indr,indc);
X_local = X(indr,indc);
Y_local = Y(indr,indc);

variant_local(boundary_local==1) = 0;

m = exx_local;
m = m(~isnan(m));
clim_exx = quantile(m,[0.003,0.997]);   % clim for exx data

limit_x_low = min(X_local(:));
limit_x_high = max(X_local(:));
limit_y_low = min(Y_local(:));
limit_y_high = max(Y_local(:));

% plot
% strain map black white as background
close all;
figure;
ax1 = axes;
surface(X_local, Y_local, exx_local,'edgecolor','none');
caxis(clim_exx);
hold on;

boundary_local(boundary_local==0) = nan;
boundary_local = boundary_local * 1;
surface(X_local, Y_local, boundary_local);
alpha(ax1, 1);
colormap(ax1, 'gray');

ax2 = axes;
ax2.Visible = 'off';

surface(X_local, Y_local, variant_local,'edgecolor','none');
alpha_mat = zeros(size(variant_local));
alpha_mat(variant_local>0) = 1;
alpha(ax2,alpha_mat);
colormap(ax2, 'parula');

axes(ax1);axis equal;
axes(ax2);axis equal;

set([ax1,ax2],'Ydir','reverse','xLim',[limit_x_low,limit_x_high],'yLim',[limit_y_low,limit_y_high]);

% linkaxes([ax1,ax2]);
ax1.XTick = [];
ax1.YTick = [];


cb1 = colorbar(ax1,'Position',[.13 .11 .05 .815]);
cb2 = colorbar(ax2,'Position',[.85 .11 .05 .815]);
% pos = get(ax1, 'position');
% set([ax1,ax2], 'position',pos);

caxis(ax2, [-0.5, 6.5]);
colormap(ax2, [0.4 0.4 0.4; parula(6)]);
set(cb2,'limits',[0.5, 6.5]);
set(cb1,'fontsize',16);
title(cb1,'\epsilon_x_x');
set(cb2,'fontsize',16);
title(cb2,'Twin Variant');

print(fullfile(output_dir, 'fig graphical abstract.tiff'), '-dtiff');