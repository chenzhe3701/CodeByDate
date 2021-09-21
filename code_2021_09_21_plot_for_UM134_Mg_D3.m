%% plot the merged DIC strain maps for sample UM134_Mg_D3

data_dir = 'E:\zhec umich Drive\2021-09-11 UM134_Mg_D3 insitu SEM-DIC\SEM Data\stitched_DIC_builtin';
output_dir = 'E:\zhec umich Drive\2021-09-11 UM134_Mg_D3 insitu SEM-DIC\Analysis\';
mkdir(output_dir);

%%
for iE = 0:7
    load(fullfile(data_dir, ['global_method_stitched_',num2str(iE),'.mat']));
    exx(sigma==-1)=nan;
    myplot(x,y,exx);
    caxis([-0.14, 0.02]);
    print(fullfile(output_dir, ['exx_iE=',num2str(iE),'.tif']),'-dtiff');
    close all;
end

%% 
make_movie_using_images('folder','E:\zhec umich Drive\2021-09-11 UM134_Mg_D3 insitu SEM-DIC\Analysis', ...
    'img_prefix','exx_iE=','img_suffix','.tif');