%% Suggest to run step by step
output_dir = 'E:\zhec umich Drive\0_temp_output\all curve yield stress';
mkdir(output_dir);

%% [Task 1] plot all the curves on one plot and compare

close all;
figure; 
set(gcf,'position',[100,100,1050,650]);
hold on;
legend_str = {};
%% Mg4Al, fine grain
p = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu curve';
f = 'Mg4Al_C1_processed_loading_data.mat';
% d = matfile(fullfile(p,f));
% stress = d.stress;
% strain = d.strain;
% plot(strain, stress, '-', 'color', 'k');
% legend_str = [legend_str, 'Mg4Al C1 fine grain'];

%% 
p = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 insitu curve\borrow some data';
f = 'Mg4Al_U2_processed_loading_data';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
plot(strain, stress, '-', 'color', 'k');
legend_str = [legend_str, 'Mg4Al U2 fine grain'];
%%
p = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu curve';
f = 'Mg4Al_C3_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
plot(strain, stress, '--', 'color', 'k');
legend_str = [legend_str, 'Mg4Al C3 fine grain'];

%% Mg4Al, a little bit larger grain
p = 'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu curve';
f = 'Mg4Al_A1_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
plot(strain, stress, '-.', 'color', 'k');
legend_str = [legend_str, 'Mg4Al A1 fine grain'];
%%
p = 'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu curve';
f = 'Mg4Al_A2_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
plot(strain, stress, ':', 'color', 'k');
legend_str = [legend_str, 'Mg4Al A2 fine grain'];

%% Mg4Al, coarse grain
p = 'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu curve';
f = 'Mg4Al_B1_processed_loading_data';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
plot(strain, stress, '-', 'color', 'r');
legend_str = [legend_str, 'Mg4Al B1 coarse grain'];
%%
p = 'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu curve';
f = 'Mg4Al_B2_processed_loading_data';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
plot(strain, stress, '--', 'color', 'r');
legend_str = [legend_str, 'Mg4Al B2 coarse grain'];

%% Pure Mg, fine grain
% p = 'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu curve';
% f = 'UM134_Mg_C1_processed_loading_data.mat';
%%
p = 'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu curve';
f = 'UM134_Mg_C2_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
plot(strain, stress, '-', 'color', 'b');
legend_str = [legend_str, 'Mg UM134 C2 fine grain'];
%%
p = 'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu curve';
f = 'UM134_Mg_C3_processed_loading_data';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
plot(strain, stress, '--', 'color', 'b');
legend_str = [legend_str, 'Mg UM134 C3 fine grain'];


%% Pure Mg, coarse grain
p = 'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu curve';
f = 'UM129_Mg_C1_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
plot(strain, stress, '-', 'color', 'm');
legend_str = [legend_str, 'Mg UM129 C1 coarse grain'];
%%
p = 'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu curve';
f = 'UM129_Mg_C2_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
ind_stop = d.ind_stop;
stress = stress(1:ind_stop(14));
strain = strain(1:ind_stop(14));
plot(strain, stress, '--', 'color', 'm');
legend_str = [legend_str, 'Mg UM129 C2 coarse grain'];
%%
p = 'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu curve';
f = 'UM129_Mg_C3_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
plot(strain, stress, ':', 'color', 'm');
legend_str = [legend_str, 'Mg UM129 C3 coarse grain'];


%%
set(gca,'fontsize',16);
xlabel('Strain'); ylabel('Stress');
legend(legend_str, 'location','eastoutside');
print(fullfile(output_dir, 'Mg all curves.tif'),'-dtiff');



%% [task 2] plot individual curves, and label the yield stress
%% Mg4Al, fine grain
% p = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu curve';
% f = 'Mg4Al_C1_processed_loading_data.mat';
%% 
p = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 insitu curve\borrow some data';
f = 'Mg4Al_U2_processed_loading_data';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
ind_stop = d.ind_stop;

close all;
figure; hold on;
plot(strain(1:ind_stop(2)), stress(1:ind_stop(2)), '-', 'color', 'k');
xlabel('Strain');
ylabel('Stress');
set(gca,'fontsize',16);
title('Mg4Al U2', 'fontweight', 'normal');

% fit 
ind = find(stress<-80,1,'first');
y = stress(1:ind);
x = strain(1:ind);
mdl = fitlm(x,y);
b = mdl.Coefficients.Estimate(1);
k = mdl.Coefficients.Estimate(2);
s = sprintf('y(MPa) = %.1f x + %.1f',k,b);

fplot(@(x) k*(x+0.002)+b, [-0.004,-0.002]);
legend(s);

%%
p = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu curve';
f = 'Mg4Al_C3_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
ind_stop = d.ind_stop;

close all;
figure; hold on;
plot(strain(1:ind_stop(2)), stress(1:ind_stop(2)), '-', 'color', 'k');
xlabel('Strain');
ylabel('Stress');
set(gca,'fontsize',16);
title('Mg4Al C3', 'fontweight', 'normal');

% fit 
ind = find(stress<-80,1,'first');
y = stress(1:ind);
x = strain(1:ind);
mdl = fitlm(x,y);
b = mdl.Coefficients.Estimate(1);
k = mdl.Coefficients.Estimate(2);
s = sprintf('y(MPa) = %.1f x + %.1f',k,b);

fplot(@(x) k*(x+0.002)+b, [-0.004,-0.002]);
legend(s);


%% Mg4Al, medium grain
p = 'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu curve';
f = 'Mg4Al_A1_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
ind_stop = d.ind_stop;

close all;
figure; hold on;
plot(strain(1:ind_stop(2)), stress(1:ind_stop(2)), '-', 'color', 'k');
xlabel('Strain');
ylabel('Stress');
set(gca,'fontsize',16);
title('Mg4Al A1', 'fontweight', 'normal');

% fit 
ind = find(stress<-80,1,'first');
y = stress(1:ind);
x = strain(1:ind);
mdl = fitlm(x,y);
b = mdl.Coefficients.Estimate(1);
k = mdl.Coefficients.Estimate(2);
s = sprintf('y(MPa) = %.1f x + %.1f',k,b);

fplot(@(x) k*(x+0.002)+b, [-0.004,-0.002]);
legend(s);
%%
p = 'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu curve';
f = 'Mg4Al_A2_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
ind_stop = d.ind_stop;

close all;
figure; hold on;
plot(strain(1:ind_stop(2)), stress(1:ind_stop(2)), '-', 'color', 'k');
xlabel('Strain');
ylabel('Stress');
set(gca,'fontsize',16);
title('Mg4Al A2', 'fontweight', 'normal');

% fit 
ind = find(stress<-80,1,'first');
y = stress(1:ind);
x = strain(1:ind);
mdl = fitlm(x,y);
b = mdl.Coefficients.Estimate(1);
k = mdl.Coefficients.Estimate(2);
s = sprintf('y(MPa) = %.1f x + %.1f',k,b);

fplot(@(x) k*(x+0.002)+b, [-0.004,-0.002]);
legend(s);


%% Mg4Al, coarse grain
p = 'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu curve';
f = 'Mg4Al_B1_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
ind_stop = d.ind_stop;

close all;
figure; hold on;
plot(strain(1:ind_stop(2)), stress(1:ind_stop(2)), '-', 'color', 'k');
xlabel('Strain');
ylabel('Stress');
set(gca,'fontsize',16);
title('Mg4Al B1', 'fontweight', 'normal');

% fit 
ind = find(stress<-60,1,'first');
y = stress(1:ind);
x = strain(1:ind);
mdl = fitlm(x,y);
b = mdl.Coefficients.Estimate(1);
k = mdl.Coefficients.Estimate(2);
s = sprintf('y(MPa) = %.1f x + %.1f',k,b);

fplot(@(x) k*(x+0.002)+b, [-0.004,-0.002]);
legend(s);
%%
p = 'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu curve';
f = 'Mg4Al_B2_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
ind_stop = d.ind_stop;

close all;
figure; hold on;
plot(strain(1:ind_stop(2)), stress(1:ind_stop(2)), '-', 'color', 'k');
xlabel('Strain');
ylabel('Stress');
set(gca,'fontsize',16);
title('Mg4Al B2', 'fontweight', 'normal');

% fit 
ind = find(stress<-55,1,'first');
y = stress(1:ind);
x = strain(1:ind);
mdl = fitlm(x,y);
b = mdl.Coefficients.Estimate(1);
k = mdl.Coefficients.Estimate(2);
s = sprintf('y(MPa) = %.1f x + %.1f',k,b);

fplot(@(x) k*(x+0.002)+b, [-0.003,-0.002]);
legend(s);



%% Pure Mg, fine grain
p = 'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu curve';
f = 'UM134_Mg_C1_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
ind_stop = d.ind_stop;

close all;
figure; hold on;
plot(strain(1:ind_stop(2)), stress(1:ind_stop(2)), '-', 'color', 'k');
xlabel('Strain');
ylabel('Stress');
set(gca,'fontsize',16);
title('UM134 Mg C1', 'fontweight', 'normal');

% fit 
ind = find(stress<-40,1,'first');
y = stress(1:ind);
x = strain(1:ind);
mdl = fitlm(x,y);
b = mdl.Coefficients.Estimate(1);
k = mdl.Coefficients.Estimate(2);
s = sprintf('y(MPa) = %.1f x + %.1f',k,b);

fplot(@(x) k*(x+0.002)+b, [-0.004,-0.002]);
legend(s);
%%
p = 'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu curve';
f = 'UM134_Mg_C2_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
ind_stop = d.ind_stop;

close all;
figure; hold on;
plot(strain(1:ind_stop(2)), stress(1:ind_stop(2)), '-', 'color', 'k');
xlabel('Strain');
ylabel('Stress');
set(gca,'fontsize',16);
title('UM134 Mg C2', 'fontweight', 'normal');

% fit 
ind = find(stress<-40,1,'first');
y = stress(1:ind);
x = strain(1:ind);
mdl = fitlm(x,y);
b = mdl.Coefficients.Estimate(1);
k = mdl.Coefficients.Estimate(2);
s = sprintf('y(MPa) = %.1f x + %.1f',k,b);

fplot(@(x) k*(x+0.002)+b, [-0.005,-0.002]);
legend(s);

%%
p = 'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu curve';
f = 'UM134_Mg_C3_processed_loading_data';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
ind_stop = d.ind_stop;

close all;
figure; hold on;
plot(strain(1:ind_stop(2)), stress(1:ind_stop(2)), '-', 'color', 'k');
xlabel('Strain');
ylabel('Stress');
set(gca,'fontsize',16);
title('UM134 Mg C3', 'fontweight', 'normal');

% fit 
ind = find(stress<-40,1,'first');
y = stress(1:ind);
x = strain(1:ind);
mdl = fitlm(x,y);
b = mdl.Coefficients.Estimate(1);
k = mdl.Coefficients.Estimate(2);
s = sprintf('y(MPa) = %.1f x + %.1f',k,b);

fplot(@(x) k*(x+0.002)+b, [-0.004,-0.002]);
legend(s);


%% Pure Mg, coarse grain
p = 'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu curve';
f = 'UM129_Mg_C1_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
ind_stop = d.ind_stop;

close all;
figure; hold on;
plot(strain(1:ind_stop(2)), stress(1:ind_stop(2)), '-', 'color', 'k');
xlabel('Strain');
ylabel('Stress');
set(gca,'fontsize',16);
title('UM129 Mg C1', 'fontweight', 'normal');

% fit 
ind = find(stress<-25,1,'first');
y = stress(1:ind);
x = strain(1:ind);
mdl = fitlm(x,y);
b = mdl.Coefficients.Estimate(1);
k = mdl.Coefficients.Estimate(2);
s = sprintf('y(MPa) = %.1f x + %.1f',k,b);

fplot(@(x) k*(x+0.002)+b, [-0.003,-0.002]);
legend(s);
%%
p = 'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu curve';
f = 'UM129_Mg_C2_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
ind_stop = d.ind_stop;
stress = stress(1:ind_stop(14));
strain = strain(1:ind_stop(14));
ind_stop = d.ind_stop;

close all;
figure; hold on;
plot(strain(1:ind_stop(2)), stress(1:ind_stop(2)), '-', 'color', 'k');
xlabel('Strain');
ylabel('Stress');
set(gca,'fontsize',16);
title('UM129 Mg C2', 'fontweight', 'normal');

% fit 
ind = find(stress<-21,1,'first');
y = stress(1:ind);
x = strain(1:ind);
mdl = fitlm(x,y);
b = mdl.Coefficients.Estimate(1);
k = mdl.Coefficients.Estimate(2);
s = sprintf('y(MPa) = %.1f x + %.1f',k,b);

fplot(@(x) k*(x+0.002)+b, [-0.003,-0.002]);
legend(s);
%%
p = 'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu curve';
f = 'UM129_Mg_C3_processed_loading_data.mat';
d = matfile(fullfile(p,f));
stress = d.stress;
strain = d.strain;
ind_stop = d.ind_stop;

close all;
figure; hold on;
plot(strain(1:ind_stop(2)), stress(1:ind_stop(2)), '-', 'color', 'k');
xlabel('Strain');
ylabel('Stress');
set(gca,'fontsize',16);
title('UM129 Mg C3', 'fontweight', 'normal');

% fit 
ind = find(stress<-22,1,'first');
y = stress(1:ind);
x = strain(1:ind);
mdl = fitlm(x,y);
b = mdl.Coefficients.Estimate(1);
k = mdl.Coefficients.Estimate(2);
s = sprintf('y(MPa) = %.1f x + %.1f',k,b);

fplot(@(x) k*(x+0.002)+b, [-0.0028,-0.002]);
legend(s);


