%% [Fig 5] twin pct vs global strain, selected samples

output_dir = 'E:\zhec umich Drive\0_temp_output\2022 gs alloy paper';
mkdir(output_dir);

%% Mg4Al samples
dir_Mg4Al_C3 = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis';
dir_Mg4Al_U2 = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis';

%% Mg4Al medium grain
dir_Mg4Al_A1 = 'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis';
dir_Mg4Al_A2 = 'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis';

%% Mg4Al coarse grain
dir_Mg4Al_B1 = 'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis';
dir_Mg4Al_B2 = 'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis';

%% UM134 unalloyed Mg samples
dir_UM134_Mg_C2 = 'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis';
dir_UM134_Mg_C3 = 'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis';

% UM129 unalloyed Mg (coarse grain) samples
dir_UM129_Mg_C2 = 'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis';
dir_UM129_Mg_C3 = 'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis';

ylim = [-2 65];
xlim = [-0.05, 0.005];

colors = mix_colorA(2);
%% For Mg4Al alloys
close all;
figure; 
hold on;

load(fullfile(dir_Mg4Al_C3, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, 'o:', 'linewidth',1.5,'markersize',8, 'color',colors(1,:), 'markerfacecolor',colors(1,:));

load(fullfile(dir_Mg4Al_U2, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, 's:', 'linewidth',1.5,'markersize',8, 'color',colors(2,:), 'markerfacecolor',colors(2,:));

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',xlim,'ylim',ylim,'fontsize',18,'fontweight','normal');
xlabel('Strain (EBSD Estimated)');
ylabel('Twin Area Percent (%)');
legend('Mg4Al C3 (fine grain)','Mg4Al U2 (fine grain)');

print(fullfile(output_dir, 'Mg4Al fine all twin pct.tiff'),'-dtiff');

%% For Mg4Al medium size grain
close all;
figure; 
hold on;

load(fullfile(dir_Mg4Al_A1, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, 'o:', 'linewidth',1.5,'markersize',8, 'color',colors(1,:), 'markerfacecolor',colors(1,:));

load(fullfile(dir_Mg4Al_A2, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, 's:', 'linewidth',1.5,'markersize',8, 'color',colors(2,:), 'markerfacecolor',colors(2,:));

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',xlim,'ylim',ylim,'fontsize',18,'fontweight','normal');
xlabel('Strain (EBSD Estimated)');
ylabel('Twin Area Percent (%)');
legend('Mg4Al A1 (medium grain)','Mg4Al A2 (medium grain)');

print(fullfile(output_dir, 'Mg4Al medium all twin pct.tiff'),'-dtiff');

%% For Mg4Al coarse grain alloys
close all;
figure; 
hold on;

load(fullfile(dir_Mg4Al_B1, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, 'o:', 'linewidth',1.5,'markersize',8, 'color',colors(1,:), 'markerfacecolor',colors(1,:));

load(fullfile(dir_Mg4Al_B2, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, 's:', 'linewidth',1.5,'markersize',8, 'color',colors(2,:), 'markerfacecolor',colors(2,:));

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',xlim,'ylim',ylim,'fontsize',18,'fontweight','normal');
xlabel('Strain (EBSD Estimated)');
ylabel('Twin Area Percent (%)');
legend('Mg4Al B1 (coarse grain)','Mg4Al B2 (coarse grain)');

print(fullfile(output_dir, 'Mg4Al coarse all twin pct.tiff'),'-dtiff');


%% For UM134 alloys
figure; 
hold on;

load(fullfile(dir_UM134_Mg_C2, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, 'o:', 'linewidth',1.5,'markersize',8, 'color',colors(1,:), 'markerfacecolor',colors(1,:));

load(fullfile(dir_UM134_Mg_C3, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, 's:', 'linewidth',1.5,'markersize',8, 'color',colors(2,:), 'markerfacecolor',colors(2,:));

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',xlim,'ylim',ylim,'fontsize',18,'fontweight','normal');
xlabel('Strain (EBSD Estimated)');
ylabel('Twin Area Percent (%)');
legend('Mg UM134 C2 (fine grain)','Mg UM134 C3 (fine grain)');

print(fullfile(output_dir, 'UM134 Mg fine all twin pct.tiff'),'-dtiff');

%% For UM129 coarse grain Mg
figure; 
hold on;

load(fullfile(dir_UM129_Mg_C2, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, 'o:', 'linewidth',1.5,'markersize',8, 'color',colors(1,:), 'markerfacecolor',colors(1,:));

load(fullfile(dir_UM129_Mg_C3, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
errorbar(strain_ebsd, tAvg, tStd, 's:', 'linewidth',1.5,'markersize',8, 'color',colors(2,:), 'markerfacecolor',colors(2,:));

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',xlim,'ylim',ylim,'fontsize',18,'fontweight','normal');
xlabel('Strain (EBSD Estimated)');
ylabel('Twin Area Percent (%)');
legend('Mg UM129 C2 (coarse grain)', 'Mg UM129 C3 (coarse grain)');

print(fullfile(output_dir, 'UM129 Mg coarse all twin pct.tiff'),'-dtiff');

%% 
close all;


