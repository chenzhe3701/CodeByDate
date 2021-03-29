
%% output dir
output_dir = 'E:\zhec umich Drive\0_temp_output';

%% Mg4Al_S1
file_name = 'E:\Mg4Al_S1_insitu\EBSD Data\Mg4Al s1.osc';
pole_figure_by_mtex(file_name);
print(fullfile(output_dir,'Mg4Al_S1 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% Mg4Al_C1
file_name = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD\Mg4Al_C1_iE=0.osc';
pole_figure_by_mtex(file_name);
print(fullfile(output_dir,'Mg4Al_C1 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% Mg4Al_C3
file_name = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\Mg4Al_C3 iE=0.osc';
pole_figure_by_mtex(file_name);
print(fullfile(output_dir,'Mg4Al_C3 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% Mg4Al_U2
file_name = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\Mg4Al_U2 iE=0.osc';
pole_figure_by_mtex(file_name);
print(fullfile(output_dir,'Mg4Al_U2 pole figure.tiff'),'-dtiff');
winopen(output_dir);


%% UM134_Mg_C2
file_name = 'E:\zhec umich Drive\2021-01-15 UM134 Mg_C2 insitu EBSD\UM134_Mg_C2 iE=0.osc';
pole_figure_by_mtex(file_name);
print(fullfile(output_dir,'UM134_Mg_C2 pole figure.tiff'),'-dtiff');
winopen(output_dir);

%% UM134_Mg_C3
file_name = 'E:\zhec umich Drive\2021-01-29 UM134 Mg_C3 insitu EBSD\UM134_Mg_C3 iE=0.osc';
pole_figure_by_mtex(file_name);
print(fullfile(output_dir,'UM134_Mg_C3 pole figure.tiff'),'-dtiff');
winopen(output_dir);

