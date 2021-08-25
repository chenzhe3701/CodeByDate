% combine available samples

output_dir = 'E:\zhec umich Drive\0_temp_output\all twin pct';
mkdir(output_dir);

%% Mg4Al samples
dir_Mg4Al_S1 = 'E:\Mg4Al_S1_insitu\Summary';
dir_Mg4Al_C1 = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD\analysis';
dir_Mg4Al_C3 = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis';
dir_Mg4Al_U2 = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis';

% UM134 unalloyed Mg samples
dir_UM134_Mg_C1 = 'E:\zhec umich Drive\2020-12-05 UM134 Mg_C1 insitu EBSD\analysis';
dir_UM134_Mg_C2 = 'E:\zhec umich Drive\2021-01-15 UM134 Mg_C2 insitu EBSD\analysis';
dir_UM134_Mg_C3 = 'E:\zhec umich Drive\2021-01-29 UM134 Mg_C3 insitu EBSD\analysis';

% UM129 unalloyed Mg (coarse grain) samples
dir_UM129_Mg_C1 = 'E:\zhec umich Drive\2021-06-29 UM129 Mg_C1 insitu EBSD\analysis';
dir_UM129_Mg_C2 = 'E:\zhec umich Drive\2021-08-20 UM129 Mg_C2 insitu EBSD\analysis';

ylim = [-2 60];
xlim = [-0.05, 0.005];

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
set(gca,'xlim',xlim,'ylim',ylim,'fontsize',18,'fontweight','normal');
xlabel('Strain from ebsd (or DIC) estimate');
ylabel('Twin Area Percent (%)');
legend('Mg4Al S1','Mg4Al C1','Mg4Al C3','Mg4Al U2');

print(fullfile(output_dir, 'Mg4Al all twin pct.tiff'),'-dtiff');

%% For UM134 alloys
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
set(gca,'xlim',xlim,'ylim',ylim,'fontsize',18,'fontweight','normal');
xlabel('Strain from ebsd estimate');
ylabel('Twin Area Percent (%)');
legend('UM134 Mg C1','UM134 Mg C2','UM134 Mg C3');

print(fullfile(output_dir, 'UM134 Mg all twin pct.tiff'),'-dtiff');

%% For UM129 coarse grain Mg
figure; 
hold on;

load(fullfile(dir_UM129_Mg_C1, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, '.-', 'linewidth',1.5,'markersize',24);

load(fullfile(dir_UM129_Mg_C2, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, '.-', 'linewidth',1.5,'markersize',24);

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',xlim,'ylim',ylim,'fontsize',18,'fontweight','normal');
xlabel('Strain from ebsd estimate');
ylabel('Twin Area Percent (%)');
legend('UM129 Mg C1 (coarse grain)', 'UM129 Mg C2 (coarse grain)');

print(fullfile(output_dir, 'UM129 Mg all twin pct.tiff'),'-dtiff');

%% 
close all;


