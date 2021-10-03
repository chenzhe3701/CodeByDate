% plot exx,exy,eyy maps for Mg4Al Mohsen sample tested on 2021-10-01
addChenFunction;

working_dir = 'E:\zhec umich Drive\2021-10-01 MgAl insitu SEM-DIC\SEM Data\stitched DIC';

exx_mean = [];

for iE = 1:3
    d = matfile(fullfile(working_dir, ['DIC_merged_',num2str(iE),'.mat']));
    x = d.x;
    y = d.y;
    exx = d.exx;
    exy = d.exy;
    eyy = d.eyy;
    sigma = d.sigma;
    
    exx(sigma==-1) = nan;
    exy(sigma==-1) = nan;
    eyy(sigma==-1) = nan;
    
    fprintf('# of points with sigma=-1: %d\n', sum(sigma(:)==-1));
    fprintf('pct of points with sigma=-1: %0.2f %% \n', 100*sum(sigma(:)==-1)/sum(numel(sigma(:))));
    
    exx_mean(iE) = nanmean(exx(:));
    
    myplot(x,y,exx);
    xlabel('x (pixels)');
    ylabel('y (pixels)');
    title_str = sprintf('\\epsilon_{xx}, load step %d, \\epsilon^{global} = %.3f', iE, exx_mean(iE));
    title(title_str, 'fontweight','normal');
    set(gca,'fontsize',14);
    print(fullfile(working_dir, ['exx load step ',num2str(iE),'.tiff']),'-dtiff');
    
    myplot(x,y,exy);
    xlabel('x (pixels)');
    ylabel('y (pixels)');
    title_str = sprintf('\\epsilon_{xy}, load step %d, \\epsilon^{global} = %.3f', iE, exx_mean(iE));
    title(title_str, 'fontweight','normal');
    set(gca,'fontsize',14);
    print(fullfile(working_dir, ['exy load step ',num2str(iE),'.tiff']),'-dtiff');
    
    myplot(x,y,eyy);
    xlabel('x (pixels)');
    ylabel('y (pixels)');
    title_str = sprintf('\\epsilon_{yy}, load step %d, \\epsilon^{global} = %.3f', iE, exx_mean(iE));
    title(title_str, 'fontweight','normal');
    set(gca,'fontsize',14);
    print(fullfile(working_dir, ['eyy load step ',num2str(iE),'.tiff']),'-dtiff');
    
    close all;
    
    
end
