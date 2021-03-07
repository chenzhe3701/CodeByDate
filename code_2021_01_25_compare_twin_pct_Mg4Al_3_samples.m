% compare twin pct of Mg4Al
clear;clc;close all;
addChenFunction;

%% load data
load('D:\p\m\CodeByDate\data\Mg4Al_twinPct_summary.mat','tbl_s1','tbl_c1','tbl_c3');

close all;

figure; hold on;
errorbar(tbl_s1.strain_corrected, 100*tbl_s1.t_avg, 100*tbl_s1.t_std, '.-', 'color', 'r', 'linewidth',1.5,'markersize',24);
errorbar(tbl_c1.strain_corrected, 100*tbl_c1.t_avg, 100*tbl_c1.t_std, '.-', 'color', 'b', 'linewidth',1.5,'markersize',24);
errorbar(tbl_c3.strain_corrected, 100*tbl_c3.t_avg, 100*tbl_c3.t_std, '.-', 'color', 'k', 'linewidth',1.5,'markersize',24);
legend({'Mg4Al SEM-DIC', 'Mg4Al EBSD sample C1', 'Mg4Al EBSD sample C3'});
xlabel('strain, estimated by DIC or EBSD');
ylabel('twin percent (%)');
set(gca, 'fontsize', 16, 'xlim', [-0.035, 0.005], 'ylim', [0 50]);

figure; hold on;
errorbar(tbl_s1.strain_corrected, 100*tbl_s1.t_avg, 100*tbl_s1.t_std, '.-', 'color', 'r', 'linewidth',1.5,'markersize',24);
errorbar(tbl_c1.strain_corrected, 100*tbl_c1.t_avg, 100*tbl_c1.t_std, '.-', 'color', 'b', 'linewidth',1.5,'markersize',24);
errorbar(tbl_c3.strain_sg, 100*tbl_c3.t_avg, 100*tbl_c3.t_std, '.-', 'color', 'k', 'linewidth',1.5,'markersize',24);
xlabel('strain');
ylabel('twin percent');
set(gca, 'fontsize', 16, 'xlim', [-0.035, 0.005], 'ylim', [0 50]);
legend({'\fontsize{12}Mg4Al SEM-DIC, estimated from strain field', ...
    '\fontsize{12}Mg4Al EBSD C1, estimated from EBSD', ...
    '\fontsize{12}Mg4Al EBSD C3, strain gage'}, ...
    'location','best');

%% Codes below: nalysis and save data
%% Mg4Al_S1, see code '2020_03_30_summarize_Mg4Al_S1.m' for further
d = load('D:\p\m\DIC_Analysis\setting_for_real_samples\Mg4Al_S1_setting.mat');
strain_sg =         [-0.0075, -0.0150, -0.0250, -0.0230, -0.0170, 0];
strain_corrected =  [-0.0012, -0.0117, -0.0186, -0.0172, -0.0124, 0];

% load('D:\p\m\DIC_Analysis\temp_results\20200325_0506_new_variant_map_Mg4Al_S1.mat','variantMapCleanedCell');
load('E:\Mg4Al_S1_insitu\Summary\twinVariantMapCell.mat','variantMapCleanedCell');
load('E:\Mg4Al_S1_insitu\Summary\Mg4Al_S1_EBSD_organized.mat','ID');

iE_start = 1;
iE_stop = 6;
for iE = iE_start:iE_stop
    variantMap = variantMapCleanedCell{iE};
    if iE==iE_start
        [nR,nC] = size(variantMap);
    end
    nr = 3;
    nc = 3;
    for ir=1:nr
       for ic = 1:nc
           subTrueTwinMap = variantMap([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic);
           subIDMap = ID([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic);
           twinPct((ir-1)*nc+ic,iE) = sum(subTrueTwinMap(:)>0)/sum(subIDMap(:)>0);
       end
    end
end
tAvg = mean(twinPct);
tStd = std(twinPct);

tbl_s1 = table(strain_sg(:), strain_corrected(:), tAvg(:), tStd(:), ...
    'VariableNames', {'strain_sg','strain_corrected', 't_avg','t_std'})

%% Mg4Al_C1
working_dir = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD';
save_dir = [working_dir, '\analysis'];
load(fullfile(save_dir,'variant_maps.mat'),'variantMap');

strain_sg =   [0, -0.0010, -0.0080, -0.0150, -0.0250, -0.0230, -0.0170, -0.003];
strain_corrected = [0, -0.0009, -0.0054, -0.0093, -0.0108, -0.0097, -0.0024, 0.0041];
stress = [0, -78, -89, -103, -127, 10, 91, 138];
twinPct = zeros(12,8);

for iE = 2:6
    load(fullfile(save_dir,['data_with_ID_overwrite_iE_',num2str(iE),'.mat']), 'ID');
    variant_map = variantMap{iE};
    if iE==2
        [nR,nC] = size(variant_map);
    end
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

%% Mg4Al_C3
working_dir = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD';
save_dir = [working_dir, '\analysis'];
load(fullfile(save_dir,'variant_maps.mat'),'variant_point_wise');

strain_corrected = [0, -0.0108, -0.0218, -0.0300, ...
    -0.0296, -0.0241, -0.0120, -0.0025, ...
    -0.0100, -0.0192, -0.0312, ...
    -0.0225, -0.0124, -0.0036];
    
% strain from image analysis, iE=0:13 + last unloaded step
strain_sg = [0, -0.0075, -0.015, -0.025, ...
    -0.023, -0.017, -0.0075, 0, ...
    -0.0075, -0.015, -0.025, ...
    -0.017, -0.0073, 0];

twinPct = zeros(12,14);     % areas x # iEs

for iE = 1:13
    load(fullfile(save_dir,['Mg4Al_C3_grain_file_iE_',num2str(iE),'.mat']), 'ID');
    
    variant_map = variant_point_wise{iE};
    if iE==1
        [nR,nC] = size(variant_map);
    end
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

twinPct
tAvg = mean(twinPct);
tStd = std(twinPct);

tbl_c3 = table(strain_sg(:), strain_corrected(:), tAvg(:), tStd(:), ...
    'VariableNames', {'strain_sg','strain_corrected', 't_avg','t_std'})

%% Save tables to data folder, no need to analyze everytime
save('D:\p\m\CodeByDate\data\Mg4Al_twinPct_summary.mat','tbl_s1','tbl_c1','tbl_c3');




