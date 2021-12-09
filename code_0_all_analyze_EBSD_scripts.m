%% [part A] all scripts analyze EBSD
open('code_2021_03_02_Mg4Al_U2_analyze_EBSD.m');    % done
winopen('E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD');
%%
open('code_2021_02_02_Mg4Al_C3_analyze_EBSD.m');    % done
winopen('E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD');
%%
open('code_2021_04_01_Mg4Al_C1_analyze_EBSD.m');    % done
winopen('E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD');


%%
open('code_2021_11_08_Mg4Al_A1_analyze_EBSD.m');    % done  
winopen('E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD');
%%
open('code_2021_11_09_Mg4Al_A2_analyze_EBSD.m');    % done
winopen('E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD');

%% There is no individual scripts for Mg4Al_B1 and Mg4Al_B2

%%
open('code_2021_02_25_UM134_Mg_C3_analyze_EBSD.m'); % done 
winopen('E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD');
%%
open('code_2021_04_02_UM134_Mg_C1_analyze_EBSD.m'); % done
winopen('E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD');
%%
open('code_2021_04_02_UM134_Mg_C2_analyze_EBSD.m'); % done
winopen('E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD');


%%
open('code_2021_07_14_UM129_Mg_C1_analyze_EBSD.m'); % done
winopen('E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD');
%%
open('code_2021_08_22_UM129_Mg_C2_analyze_EBSD.m'); % done
winopen('E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD');
%%
open('code_2021_09_07_UM129_Mg_C3_analyze_EBSD.m'); % done
winopen('E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD');



%% [part B] all scripts analyze curve

%% Mg4Al, fine grain
p = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu curve';
f = 'Mg4Al_C1_processed_loading_data.mat';
open('code_2020_10_25_Mg4Al_C1_curve.m');
winopen(p);
%% 
p = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 insitu curve';
f = 'Mg4Al_U2_processed_loading_data';
open('code_2021_03_01_Mg4Al_U2_curve.m');
winopen(p);
%%
p = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu curve';
f = 'Mg4Al_C3_processed_loading_data.mat';
open('code_2020_12_29_Mg4Al_C3_curve.m');
winopen(p);


%% Mg4Al, medium grain
p = 'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu curve';
f = 'Mg4Al_A1_processed_loading_data.mat';
open('code_2021_11_08_Mg4Al_A1_curve.m');
winopen(p);
%%
p = 'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu curve';
f = 'Mg4Al_A2_processed_loading_data.mat';
open('code_2021_11_08_Mg4Al_A2_curve.m');
winopen(p);


%% Mg4Al, coarse grain
p = 'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu curve';
f = 'Mg4Al_B1_processed_loading_data.mat';
open('code_2021_12_04_Mg4Al_B1_curve.m');
winopen(p);
%%
p = 'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu curve';
f = 'Mg4Al_B2_processed_loading_data.mat';
open('code_2021_12_06_Mg4Al_B2_curve.m');
winopen(p);


%% Pure Mg, fine grain
p = 'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu curve';
f = 'UM134_Mg_C1_processed_loading_data.mat';
open('code_2020_12_08_UM134_Mg_C1_curve.m');
winopen(p);
%%
p = 'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu curve';
f = 'UM134_Mg_C2_processed_loading_data.mat';
open('code_2021_01_17_UM134_Mg_C2_curve.m');
winopen(p);
%%
p = 'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu curve';
f = 'UM134_Mg_C3_processed_loading_data';
open('code_2021_02_02_UM134_Mg_C3_curve.m');
winopen(p);


%% Pure Mg, coarse grain
p = 'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu curve';
f = 'UM129_Mg_C1_processed_loading_data.mat';
open('code_2021_07_14_UM129_Mg_C1_curve.m');
winopen(p);
%%
p = 'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu curve';
f = 'UM129_Mg_C2_processed_loading_data.mat';
open('code_2021_08_22_UM129_Mg_C2_curve.m');
winopen(p);
%%
p = 'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu curve';
f = 'UM129_Mg_C3_processed_loading_data.mat';
open('code_2021_09_07_UM129_Mg_C3_curve.m');
winopen(p);



