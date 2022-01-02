
%% output dir
output_dir = 'E:\zhec umich Drive\0_temp_output\PF individual scale';
clims = [];

% comment to choose clims
% output_dir = 'E:\zhec umich Drive\0_temp_output\PF same scale';
% clims = {[0,7.2], [0,4.2]};

mkdir(output_dir);
mkdir(fullfile(output_dir, 'big AOI'));
mkdir(fullfile(output_dir, 'avg'));
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


%% Mg4Al_A1 (medium grain)
file_name = 'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\Mg4Al_A1 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'Mg4Al_A1 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% Mg4Al_A2 (medium grain)
file_name = 'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\Mg4Al_A2 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'Mg4Al_A2 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% Mg4Al_B1 (coarse grain)
file_name = 'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\Mg4Al_B1 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'Mg4Al_B1 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% Mg4Al_B2 (coarse grain)
file_name = 'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\Mg4Al_b2 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'Mg4Al_B2 pole figure.tiff'),'-dtiff');
winopen(output_dir);


%% UM134_Mg_C1
file_name = 'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\UM134_Mg_C1 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'UM134_Mg_C1 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% UM134_Mg_C2
file_name = 'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\UM134_Mg_C2 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'UM134_Mg_C2 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% UM134_Mg_C3
file_name = 'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\UM134_Mg_C3 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'UM134_Mg_C3 pole figure.tiff'),'-dtiff');
winopen(output_dir);


%% UM129_Mg_C1
file_name = 'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\UM129_Mg_C1 EBSD iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'UM129_Mg_C1 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% UM129_Mg_C2
file_name = 'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\UM129_Mg_C2 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'UM129_Mg_C2 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% UM129_Mg_C3
file_name = 'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\UM129_Mg_C3 iE=0.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir,'UM129_Mg_C3 pole figure.tiff'),'-dtiff');
winopen(output_dir);


%% big aoi
file_name_Mg4Al_C2_big = {'E:\zhec umich Drive\2020-11-12 Mg4Al_C2 EBSD big AOI\Mg4Al_C2 big aoi left.osc', ...
	'E:\zhec umich Drive\2020-11-12 Mg4Al_C2 EBSD big AOI\Mg4Al_C2 big aoi mid.osc', ...
    'E:\zhec umich Drive\2020-11-12 Mg4Al_C2 EBSD big AOI\Mg4Al_C2 big aoi right.osc'};
pole_figure_by_mtex_multi_file(file_name_Mg4Al_C2_big, 'setting_str','setting 1', 'clims',clims);
print(fullfile(output_dir, 'big AOI', 'Mg4Al_C2 big aois avg PF.tiff'), '-dtiff');

file_name_UM134_Mg_C1_big = {'E:\zhec umich Drive\2020-12-01 UM134_Mg_C1 big AOI\UM134_Mg_C1 big aoi left.osc', ...
    'E:\zhec umich Drive\2020-12-01 UM134_Mg_C1 big AOI\UM134_Mg_C1 big aoi mid.osc', ...
    'E:\zhec umich Drive\2020-12-01 UM134_Mg_C1 big AOI\UM134_Mg_C1 big aoi right.osc'};
pole_figure_by_mtex_multi_file(file_name_UM134_Mg_C1_big, 'setting_str','setting 1', 'clims',clims);
print(fullfile(output_dir, 'big AOI', 'UM134_Mg_C1 big aois avg PF.tiff'), '-dtiff');

file_name = 'E:\zhec umich Drive\2021-10-26 Mg4Al_A1 coarse grain big aoi EBSD\Mg4Al_A1 coarse grain big aoi.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir, 'big AOI', 'Mg4Al_A1 big aoi PF.tiff'),'-dtiff');

file_name = 'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 big aoi\Mg4Al_B1 EBSD big aoi.osc';
pole_figure_by_mtex(file_name, 'clims',clims);
print(fullfile(output_dir, 'big AOI', 'Mg4Al_B1 big aoi PF.tiff'),'-dtiff');

winopen(output_dir);

%%
close all;


%% Plot avg pole figures for different materials 
file_name_Mg4Al = {% 'E:\Mg4Al_S1_insitu\EBSD Data\Mg4Al s1.osc', ...
    'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD\Mg4Al_C1_iE=0.osc', ...
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\Mg4Al_C3 iE=0.osc', ...
    'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\Mg4Al_U2 iE=0.osc'};

file_name_Mg4Al_medium = {'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\Mg4Al_A1 iE=0.osc', ...
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\Mg4Al_A2 iE=0.osc'};

file_name_Mg4Al_coarse = {'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\Mg4Al_B1 iE=0.osc', ...
    'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\Mg4Al_B2 iE=0.osc'};

file_name_UM134_Mg = {'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\UM134_Mg_C1 iE=0.osc', ...
    'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\UM134_Mg_C2 iE=0.osc', ...
    'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\UM134_Mg_C3 iE=0.osc'};

file_name_UM129_Mg = {'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\UM129_Mg_C1 EBSD iE=0.osc', ...
    'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\UM129_Mg_C2 iE=0.osc', ...
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\UM129_Mg_C3 iE=0.osc'};

pole_figure_by_mtex_multi_file(file_name_Mg4Al, 'setting_str','setting 1', 'clims',clims);
print(fullfile(output_dir, 'avg', 'Mg4Al_fine avg pole figure.tiff'), '-dtiff');

pole_figure_by_mtex_multi_file(file_name_Mg4Al_medium, 'setting_str','setting 1', 'clims',clims);
print(fullfile(output_dir, 'avg', 'Mg4Al_medium avg pole figure.tiff'), '-dtiff');

pole_figure_by_mtex_multi_file(file_name_Mg4Al_coarse, 'setting_str','setting 1', 'clims',clims);
print(fullfile(output_dir, 'avg', 'Mg4Al_coarse avg pole figure.tiff'), '-dtiff');

pole_figure_by_mtex_multi_file(file_name_UM134_Mg, 'setting_str','setting 1', 'clims',clims);
print(fullfile(output_dir, 'avg', 'UM134_Mg_fine avg pole figure.tiff'), '-dtiff');

pole_figure_by_mtex_multi_file(file_name_UM129_Mg, 'setting_str','setting 1', 'clims',clims);
print(fullfile(output_dir, 'avg', 'UM129_Mg_coarse avg pole figure.tiff'), '-dtiff');

close all;
