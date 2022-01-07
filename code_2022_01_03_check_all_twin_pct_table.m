% summarize twin percent in a excel file

output_dir = 'E:\zhec umich Drive\0_temp_output\2022 gs alloy paper';
mkdir(output_dir);

% Mg4Al samples
dir_Mg4Al_S1 = 'E:\Mg4Al_S1_insitu\Summary';
dir_Mg4Al_C1 = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD\analysis';
dir_Mg4Al_C3 = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis';
dir_Mg4Al_U2 = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis';

% Mg4Al medium grain
dir_Mg4Al_A1 = 'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis';
dir_Mg4Al_A2 = 'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis';

% Mg4Al coarse grain
dir_Mg4Al_B1 = 'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis';
dir_Mg4Al_B2 = 'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis';

% UM134 unalloyed Mg samples
dir_UM134_Mg_C1 = 'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis';
dir_UM134_Mg_C2 = 'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis';
dir_UM134_Mg_C3 = 'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis';

% UM129 unalloyed Mg (coarse grain) samples
dir_UM129_Mg_C1 = 'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis';
dir_UM129_Mg_C2 = 'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis';
dir_UM129_Mg_C3 = 'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis';

ylim = [-2 65];
xlim = [-0.05, 0.005];

output_file = fullfile(output_dir, 'test.xlsx');

%% For Mg4Al alloys

load(fullfile(dir_Mg4Al_C3, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
t = table(tbl.iE, tAvg, 'VariableNames', {'iE', 'Mg4Al_C3'});
T = t;

load(fullfile(dir_Mg4Al_U2, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
t = table(tbl.iE, tAvg, 'VariableNames', {'iE', 'Mg4Al_U2'});
T = innerjoin(T, t, 'Keys','iE');

% For Mg4Al medium size grain
load(fullfile(dir_Mg4Al_A1, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
t = table(tbl.iE, tAvg, 'VariableNames', {'iE', 'Mg4Al_A1'});
T = innerjoin(T, t, 'Keys','iE');

load(fullfile(dir_Mg4Al_A2, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
t = table(tbl.iE, tAvg, 'VariableNames', {'iE', 'Mg4Al_A2'});
T = innerjoin(T, t, 'Keys','iE');


% For Mg4Al coarse grain alloys
load(fullfile(dir_Mg4Al_B1, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
t = table(tbl.iE, tAvg, 'VariableNames', {'iE', 'Mg4Al_B1'});
T = innerjoin(T, t, 'Keys','iE');

load(fullfile(dir_Mg4Al_B2, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
t = table(tbl.iE, tAvg, 'VariableNames', {'iE', 'Mg4Al_B2'});
T = innerjoin(T, t, 'Keys','iE');


% For UM134 alloys
load(fullfile(dir_UM134_Mg_C2, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
t = table(tbl.iE, tAvg, 'VariableNames', {'iE', 'UM134_Mg_C2'});
T = innerjoin(T, t, 'Keys','iE');

load(fullfile(dir_UM134_Mg_C3, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
t = table(tbl.iE, tAvg, 'VariableNames', {'iE', 'UM134_Mg_C3'});
T = innerjoin(T, t, 'Keys','iE');


% For UM129 coarse grain Mg
load(fullfile(dir_UM129_Mg_C2, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
t = table(tbl.iE, tAvg, 'VariableNames', {'iE', 'UM129_Mg_C2'});
T = innerjoin(T, t, 'Keys','iE');

load(fullfile(dir_UM129_Mg_C3, 'twin_pct'), 'tbl');
strain_ebsd = tbl.strain_ebsd;
tAvg = tbl.("twinPct %");
tStd = tbl.("twinStd %");
t = table(tbl.iE, tAvg, 'VariableNames', {'iE', 'UM129_Mg_C3'});
T = innerjoin(T, t, 'Keys','iE');


writetable(T, output_file);

