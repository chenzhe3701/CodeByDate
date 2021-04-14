% combine available samples

% Mg4Al samples
dir_Mg4Al_S1 = 'E:\Mg4Al_S1_insitu\Summary';
dir_Mg4Al_C1 = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD\analysis';
dir_Mg4Al_C3 = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis';
dir_Mg4Al_U2 = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis';

% UM134 unalloyed Mg samples
dir_UM134_Mg_C1 = 'E:\zhec umich Drive\2020-12-05 UM134 Mg_C1 insitu EBSD\analysis';
dir_UM134_Mg_C2 = 'E:\zhec umich Drive\2021-01-15 UM134 Mg_C2 insitu EBSD\analysis';
dir_UM134_Mg_C3 = 'E:\zhec umich Drive\2021-01-29 UM134 Mg_C3 insitu EBSD\analysis';

%% For Mg4Al alloys
close all;
figure; 
hold on;

load(fullfile(dir_Mg4Al_S1, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_dic;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, '.-', 'linewidth',1.5,'markersize',24);

load(fullfile(dir_Mg4Al_C1, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, '.-', 'linewidth',1.5,'markersize',24);

load(fullfile(dir_Mg4Al_C3, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, '.-', 'linewidth',1.5,'markersize',24);

load(fullfile(dir_Mg4Al_U2, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, '.-', 'linewidth',1.5,'markersize',24);

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',[-0.05, 0.005],'ylim',[-2 50],'fontsize',18,'fontweight','normal');
xlabel('Strain from ebsd (or DIC) estimate');
ylabel('Twin Area Percent (%)');
legend('Mg4Al S1','Mg4Al C1','Mg4Al C3','Mg4Al U2');

%% For UM134 alloys
close all;
figure; 
hold on;

load(fullfile(dir_UM134_Mg_C1, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, '.-', 'linewidth',1.5,'markersize',24);

load(fullfile(dir_UM134_Mg_C2, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, '.-', 'linewidth',1.5,'markersize',24);

load(fullfile(dir_UM134_Mg_C3, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, '.-', 'linewidth',1.5,'markersize',24);

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',[-0.05, 0.005],'ylim',[-2 35],'fontsize',18,'fontweight','normal');
xlabel('Strain from ebsd estimate');
ylabel('Twin Area Percent (%)');
legend('UM134 Mg C1','UM134 Mg C2','UM134 Mg C3');

