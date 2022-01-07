% check all sample twin percent, without error bar
clear;clc;
addChenFunction;
% [data_dir, sample_name, sample_ID, plot_symbol, group_number]
cells = {'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1;
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1;
    'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2;
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2;
    'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3;
    'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3;
    'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'o', 4;
    'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 's', 4;
    'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 'x', 4;
    'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'o', 5;
    'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 's', 5;
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 'x', 5};

%%
clear strain twinPct
% iR = iE, iC = iSample
for ic = 1:length(cells)
    data_dir = cells{ic,1};
    sample_name = cells{ic,2};
    
    load(fullfile(data_dir,'variant_maps.mat'),'variant_point_wise');
    load(fullfile(data_dir, 'twin_pct'), 'tbl');
    strain(:,ic) = tbl.strain_ebsd;
    
    for iE = 0:13
        iB = iE + 1;
        
        load(fullfile(data_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), 'ID');
        variant_map = variant_point_wise{iB};                 
                
        twinPct(iB,ic) = sum(variant_map(:)>0)/sum(ID(:)>0)    
    end
end

%% plot twin pct vs. iE
close all;
colors = [0 0 0; mix_color(4)];

figure; hold on;
for ii = [1:6,8,9,11,12]
    symble_str = [cells{ii,4},':'];
    ig = cells{ii,5};
    plot(0:13, twinPct(:,ii), symble_str, 'color',colors(ig,:), 'linewidth',1.5, 'markersize',6, 'markerfacecolor',colors(ig,:));
end
xlabel('iE');
ylabel('Twin Fraction');

%% plot twin pct vs. strian_EBSD
close all;
colors = [0 0 0; mix_color(4)];

figure; hold on;
for ii = [1:6,8,9,11,12]
    symble_str = [cells{ii,4},':'];
    ig = cells{ii,5};
    plot(strain(:,ii), twinPct(:,ii), symble_str, 'color',colors(ig,:), 'linewidth',1.5, 'markersize',6, 'markerfacecolor',colors(ig,:));
end
xlabel('Strain (EBSD Estimate)');
ylabel('Twin Fraction');

