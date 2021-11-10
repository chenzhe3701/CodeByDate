
%% output dir
output_dir = 'E:\zhec umich Drive\0_temp_output\individual scale PF';
clims = [];

% comment to choose clims
output_dir = 'E:\zhec umich Drive\0_temp_output\same scale PF';
clims = {[0,7.2], [0,4.2]};

mkdir(output_dir);

%% Mg4Al_S1
file_name = 'E:\Mg4Al_S1_insitu\EBSD Data\Mg4Al s1.osc';
pole_figure_by_mtex(file_name, 'setting_str', 'setting 2', 'clims',clims);
print(fullfile(output_dir,'Mg4Al_S1 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% Mg4Al_C1
file_name = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD\Mg4Al_C1_iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'Mg4Al_C1 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% Mg4Al_C3
file_name = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\Mg4Al_C3 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'Mg4Al_C3 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% Mg4Al_U2
file_name = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\Mg4Al_U2 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'Mg4Al_U2 pole figure.tiff'),'-dtiff');
winopen(output_dir);


%% Mg4Al_A1 (coarse grain)
file_name = 'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\Mg4al_A1 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'Mg4Al_A1 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% Mg4Al_A2 (coarse grain)
file_name = 'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\Mg4al_A2 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'Mg4Al_A2 pole figure.tiff'),'-dtiff');
winopen(output_dir);


%% UM134_Mg_C1
file_name = 'E:\zhec umich Drive\2020-12-05 UM134 Mg_C1 insitu EBSD\UM_134 Mg_C1 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'UM134_Mg_C1 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% UM134_Mg_C2
file_name = 'E:\zhec umich Drive\2021-01-15 UM134 Mg_C2 insitu EBSD\UM134_Mg_C2 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'UM134_Mg_C2 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% UM134_Mg_C3
file_name = 'E:\zhec umich Drive\2021-01-29 UM134 Mg_C3 insitu EBSD\UM134_Mg_C3 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'UM134_Mg_C3 pole figure.tiff'),'-dtiff');
winopen(output_dir);


%% UM129_Mg_C1
file_name = 'E:\zhec umich Drive\2021-06-29 UM129 Mg_C1 insitu EBSD\UM129 Mg_C1 EBSD iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'UM129_Mg_C1 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% UM129_Mg_C2
file_name = 'E:\zhec umich Drive\2021-08-20 UM129 Mg_C2 insitu EBSD\UM129_Mg_C2 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'UM129_Mg_C2 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% UM129_Mg_C3
file_name = 'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\UM129_Mg_C3 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'UM129_Mg_C3 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%%
close all;


%% Plot avg pole figures for different materials 
file_name_Mg4Al = {% 'E:\Mg4Al_S1_insitu\EBSD Data\Mg4Al s1.osc', ...
    'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD\Mg4Al_C1_iE=0.osc', ...
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\Mg4Al_C3 iE=0.osc', ...
    'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\Mg4Al_U2 iE=0.osc'};

file_name_Mg4Al_coarse = {'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\Mg4al_A1 iE=0.osc', ...
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\Mg4al_A2 iE=0.osc'};

file_name_UM134_Mg = {'E:\zhec umich Drive\2020-12-05 UM134 Mg_C1 insitu EBSD\UM_134 Mg_C1 iE=0.osc', ...
    'E:\zhec umich Drive\2021-01-15 UM134 Mg_C2 insitu EBSD\UM134_Mg_C2 iE=0.osc', ...
    'E:\zhec umich Drive\2021-01-29 UM134 Mg_C3 insitu EBSD\UM134_Mg_C3 iE=0.osc'};

file_name_UM129_Mg = {'E:\zhec umich Drive\2021-06-29 UM129 Mg_C1 insitu EBSD\UM129 Mg_C1 EBSD iE=0.osc', ...
    'E:\zhec umich Drive\2021-08-20 UM129 Mg_C2 insitu EBSD\UM129_Mg_C2 iE=0.osc', ...
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\UM129_Mg_C3 iE=0.osc'};

pole_figure_by_mtex_multi_file(file_name_Mg4Al, 'setting_str','setting 1', 'clims',clims);
print(fullfile(output_dir, 'Mg4Al avg pole figure.tiff'), '-dtiff');

pole_figure_by_mtex_multi_file(file_name_Mg4Al_coarse, 'setting_str','setting 1', 'clims',clims);
print(fullfile(output_dir, 'Mg4Al_coarse avg pole figure.tiff'), '-dtiff');

pole_figure_by_mtex_multi_file(file_name_UM134_Mg, 'setting_str','setting 1', 'clims',clims);
print(fullfile(output_dir, 'UM134_Mg avg pole figure.tiff'), '-dtiff');

pole_figure_by_mtex_multi_file(file_name_UM129_Mg, 'setting_str','setting 1', 'clims',clims);
print(fullfile(output_dir, 'UM129_Mg avg pole figure.tiff'), '-dtiff');

close all;
