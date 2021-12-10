
clear; clc; close all;
addChenFunction;

%% load data
load('D:\p\m\CodeByDate\data\UM134_Mg_twinPct_summary.mat','tbl_c1','tbl_c2');

close all;

figure; hold on;
errorbar(tbl_c1.strain_corrected, 100*tbl_c1.t_avg, 100*tbl_c1.t_std, '.-', 'color', 'r', 'linewidth',1.5,'markersize',24);
errorbar(tbl_c2.strain_corrected, 100*tbl_c2.t_avg, 100*tbl_c2.t_std, '.-', 'color', 'b', 'linewidth',1.5,'markersize',24);
legend({'UM134 Mg sample C1', 'UM134 Mg sample C2'});
xlabel('strain (%), estimated by DIC or EBSD');
ylabel('twin percent (%)');
set(gca, 'fontsize', 16, 'xlim', [-0.032, 0.005], 'ylim', [0 30]);

figure; hold on;
errorbar(tbl_c1.strain_corrected, 100*tbl_c1.t_avg, 100*tbl_c1.t_std, '.-', 'color', 'r', 'linewidth',1.5,'markersize',24);
errorbar(tbl_c2.strain_sg, 100*tbl_c2.t_avg, 100*tbl_c2.t_std, '.-', 'color', 'b', 'linewidth',1.5,'markersize',24);
xlabel('strain (%)');
ylabel('twin percent (%)');
set(gca, 'fontsize', 16, 'xlim', [-0.032, 0.005], 'ylim', [0 30]);
legend({'\fontsize{12}UM134 Mg C1, estimated from EBSD', ...
    '\fontsize{12}UM134 Mg C2, strain gage'}, ...
    'location','best');


%% UM134_Mg_C1
working_dir = 'E:\zhec umich Drive\2020-12-05 UM_134 Mg_C1 insitu EBSD';
save_dir = [working_dir, '\analysis'];

load(fullfile(save_dir,'variant_maps.mat'),'variantMap');
% strain from strain gage, iE=0:13 + last unloaded step
strain_sg = [0, -0.0035, -0.0037, ...
    -0.0032, -0.0020, -0.0010, 0.0000, ...
    -0.0010, -0.0020, -0.0033, -0.0035, ...
    -0.0020, -0.0010, -0.0000];

% 13 load steps
strain_corrected = [1, 0.9798, 0.9759, ...
    0.9751, 0.9831, 0.9950, 1.0042, ...
    1.0001, 0.9976, 0.9926, 0.9934, ...
    0.9961, 1.0012, 1.0059] - 1;

strain_corrected = [0, -0.0202, -0.0241, ...
    -0.0249, -0.0169, -0.0050, 0.0042, ...
    0.0001, -0.0024, -0.0074, -0.0066, ...
    -0.0039, 0.0012, 0.0059];   

twinPct = zeros(14,8);

for iE = 1:13
    variant_map = variantMap{iE};
    if iE==2
        [nR,nC] = size(variant_map);
    end
    
    load(fullfile(save_dir,['data_with_ID_overwrite_iE_',num2str(iE),'.mat']), 'ID');
    
    nr = 3;
    nc = 4;
    for ir=1:nr
       for ic = 1:nc
           sub_variant_map = variant_map([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic);
           sub_ID_map = ID([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic);
           
           ii = iE + 1;
           twinPct((ir-1)*nc+ic,ii) = sum(sub_variant_map(:)>0)/sum(sub_ID_map(:)>0);
       end
    end
    
end

tAvg = mean(twinPct);
tStd = std(twinPct);

tbl_c1 = table(strain_sg(:), strain_corrected(:), tAvg(:), tStd(:), ...
    'VariableNames', {'strain_sg','strain_corrected', 't_avg','t_std'})


%% UM134_Mg_C2
working_dir = 'E:\zhec umich Drive\2021-01-15 UM_134 Mg_C2 insitu EBSD';
save_dir = [working_dir, '\analysis'];
load(fullfile(save_dir,'variant_maps.mat'),'variantMap');

% raw:
% strain_corrected = [1, 0.9844, 0.9744, 0.9681, ...
%     0.9710, 0.9783, 0.9893, 0.9978, ...
%     0.9875, 0.9788, 0.9688, ...
%     0.9796, 0.9899, 0.9986] - 1;
    
strain_corrected = [0, -0.0156, -0.0256, -0.0319, ...
    -0.0290, -0.0217, -0.0107, -0.0022, ...
    -0.0125, -0.0212, -0.0312, ...
    -0.0204, -0.0101, -0.0014];
    
% strain from strain gage, iE=0:13 + last unloaded step
strain_sg = [0, -0.0075, -0.015, -0.025, ...
    -0.023, -0.017, -0.0075, 0, ...
    -0.0075, -0.015, -0.025, ...
    -0.017, -0.0073, 0];

twinPct = zeros(12,14);

for iE = 1:13
    variant_map = variantMap{iE};
    if iE==1
        [nR,nC] = size(variant_map);
    end
   
    load(fullfile(save_dir,['UM134_Mg_C2_grain_file_iE_',num2str(iE),'.mat']), 'ID');
    
    nr = 3;
    nc = 4;    
    for ir=1:nr
       for ic = 1:nc
           sub_variant_map = variant_map([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic);
           sub_ID_map = ID([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic); % just for counting number of pixels
           
           ii = iE + 1;
           twinPct((ir-1)*nc+ic,ii) = sum(sub_variant_map(:)>0)/sum(sub_ID_map(:)>0);
       end
    end
end

tAvg = mean(twinPct);
tStd = std(twinPct);

tbl_c2 = table(strain_sg(:), strain_corrected(:), tAvg(:), tStd(:), ...
    'VariableNames', {'strain_sg','strain_corrected', 't_avg','t_std'})

%% save
save('D:\p\m\CodeByDate\data\UM134_Mg_twinPct_summary.mat','tbl_c1','tbl_c2');





