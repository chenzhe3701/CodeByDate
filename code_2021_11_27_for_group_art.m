% make some nice images: 4 load steps
% reference function: 'plot_overlay_maps_SEM_DIC'

output_dir = 'E:\zhec umich Drive\0_temp_output\mat art small';
mkdir(output_dir);

addChenFunction;

colorbarOnOff = 'on';
fontSize = 14;
c_min = -0.06;
c_max = 0.10;

%%
close all;

for iE = 7:10
    % read image
    img = imread(['E:\Ti7Al_B6_insitu_tension\SEM Data\byFov_IMAGE_DELETED_TO_REDUCE_SIZE\r2c0\Ti7Al_B6_e',num2str(iE),'_r2c0.tif']);
    % convert to index
    cmap = linspace(0,1,65536)'*[1 1 1];
    img = ind2rgb(img, cmap);
    
    % get xlim, ylim
    [x_SEM,y_SEM] = meshgrid(0:size(img,1)-1, 0:size(img,2)-1);
    x_low = min(unique(x_SEM(:)));
    x_high = max(unique(x_SEM(:)));
    y_low = min(unique(y_SEM(:)));
    y_high = max(unique(y_SEM(:)));
    
    % get strain data
    d = matfile(['E:\Ti7Al_B6_insitu_tension\SEM Data\byFov_IMAGE_DELETED_TO_REDUCE_SIZE\r2c0\Ti7Al_B6_e',num2str(iE),'_r2c0.mat']);
    x = d.x;
    y = d.y;
    u = d.u;
    v = d.v;
    sigma = d.sigma;
    exx = d.exx;
    exy = d.exy;
    eyy = d.eyy;
    xp = x+u;
    yp = y+v;
    
    val = exx;
    val(sigma==-1) = nan;

    %%    
    figure();hold on; axis off;
    
    image([x_low,x_high],[y_low,y_high],img,'AlphaData', .8);
    surf(xp,yp,val+1, val, 'FaceAlpha', .6, 'LineStyle', 'none');
    colormap('viridis');
    axis equal;
    set(gca,'xlim',[x_low,x_high],'ylim',[y_low,y_high],'ydir','reverse','xTickLabel','','yTickLabel','','xTick','','yTick','');
    view([0,90]);
    if((~exist('c_min','var'))||isempty(c_min))
        c_min = nanmean(val(:)) - 2 * nanstd(val(:));
        c_max = nanmean(val(:)) + 2 * nanstd(val(:));
    end
    caxis([c_min,c_max]);
    caxis([-0.02 0.08]);
    cbar = colorbar;
    set(cbar,'fontsize',fontSize,'visible',colorbarOnOff);
    
    xlim = [500,3500];
    ylim = [500,3500];
    set(gca,'xlim',xlim,'ylim',ylim);
    
    print(fullfile(output_dir, ['Ti7Al_B6_exx_iE=',num2str(iE),'.tif']),'-dtiff','-r600');
end

close all;

%%
%
% figure();
% ax1 = axes;
% imagesc([x_low,x_high],[y_low,y_high],img,'AlphaData', .8);
% colormap(ax1,'gray');
% axis equal;
%
% val = exx;
% val(sigma==-1)=nan;
% zp = val + 1;
% cp = val;
%
% ax2 = axes;
% ax2.Visible = 'off';
%
% % note: need to use 'facealpha'.  'alpha' as matrix does not work well.
% surface(xp,yp,zp,cp, 'edgecolor', 'none','facealpha',0.6);
% axis equal;
% set(ax1,'Ydir','reverse','xLim',[x_low,x_high],'yLim',[y_low,y_high]);
% set(ax2,'Ydir','reverse','xLim',[x_low,x_high],'yLim',[y_low,y_high]);
%
% if((~exist('c_min','var'))||isempty(c_min))
%     c_min = nanmean(exx(:)) - 2 * nanstd(exx(:));
%     c_max = nanmean(exx(:)) + 2 * nanstd(exx(:));
% end
% caxis([c_min,c_max]);
% colormap(ax2, 'parula');
%
% cbar = colorbar(ax1);
% set(cbar,'fontsize',fontSize,'visible','off');
% cbar = colorbar(ax2);
% set(cbar,'fontsize',fontSize,'visible',colorbarOnOff);
%
