%% based on script in the data folder
clear;clc;
close all;

working_dir = 'E:\Mg4Al_S1_insitu\CompressionTension Data';
output_dir = 'E:\zhec umich Drive\0_temp_output\Mg4Al_S1 analysis\loading curve';
mkdir(output_dir);

%% strain gage data
fid = fopen(fullfile(working_dir, 'Mg4Al_S1_strain_gage.TXT'));
a = textscan(fid,'%{hh:mm:ss}D %f','delimiter','\t','headerlines',6); % determined that it started 9sec earlier
ue = a{1,2};
t_sg = (0:length(ue)-1)';
fclose(fid);


%% DDS data
fid = fopen(fullfile(working_dir, 'Mg4Al_S1_insitu.TXT'));
b = textscan(fid,'%s %f %f %f %f %f','delimiter','\t','headerlines',2);
Load = b{1,3};
elongation = b{1,4};
t_dds = (0:length(Load)-1)';
fclose(fid);

A = 9.04;
stress = Load/A;

%% use for plot raw data
figure;
plot(t_sg,ue/min(ue(:)),'r');
hold;
plot(t_dds,elongation/min(elongation(:)),'b');
legend({'ue','elongation'});

%% make them the same length (device time runs in different speed ...). Note: both aligned to [dds time].
ind_ue = (6:237921)';
ind_dds = (2:237916)';

t_sg = t_sg(ind_ue,:);
ue = ue(ind_ue,:);

t_dds = t_dds(ind_dds,:);
Load = Load(ind_dds,:);
stress = stress(ind_dds,:);
elongation = elongation(ind_dds,:);

ue = round(interp1(linspace(t_dds(1),t_dds(end),length(t_sg)), ue, t_dds));
%% plot data after adjusting time
figure; hold on;
plot(t_dds,ue/min(ue(:)),'r');
plot(t_dds,elongation/min(elongation(:)),'b');
plot(t_dds,Load/min(Load(:)),'g');
legend({'ue','elongation','load'});

figure; plot(elongation,ue);
set(gca,'ylim',[-3e4, 0.5e4]);
xlabel('elongation,um');
ylabel('ue');
%% select data for active loading part
mat = [ue(:),elongation(:),stress(:)];
open mat;
%%
inds = [{1:173}, {69766:69868},  {111438:111581},...
    {143183:143233}, {174633:174860}, {206252:206543}, {237818:237900}];
rows = [];
for ii = 1:length(inds)-1
    rows = [rows,inds{ii}];
    ind_stop(ii) = length(rows);
    stress_stop(ii) = stress(rows(end));    % stress at pauses
end

% rs = [13:363, 29805:29854,  56796:56827,57292:57306,58144:58168,...
%     84932:85074, 112071:112243, 136710:137081, 161593:162218];

% reduced
time_r = t_dds(rows,:);
Load_r = Load(rows,:);
stress_r = stress(rows,:);
elongation_r = elongation(rows,:);
ue_r = ue(rows,:);

figure; plot(elongation_r, ue_r);
set(gca,'ylim',[-3e4, 0.5e4]);
xlabel('elongation,um');
ylabel('ue');
display(stress_stop);

%% No need to interp/extrap ue with elongation
% inds_interp = 260:690;
% model = fitlm(elongation_r(inds_interp,:),ue_r(inds_interp,:));
% slope = model.Coefficients.Estimate(2);
% ue_r(inds_interp(end)+1:end) = elongation_r(inds_interp(end)+1:end) * slope;
strain_r = ue_r/1000000;
%% [plot 1] stress vs gage_strain
close;
figure; hold on

plot(strain_r,stress_r,'linewidth',2);
for ii = 1:length(ind_stop)
   plot(strain_r(ind_stop(ii)),stress_r(ind_stop(ii)),'.r','markersize',18);
   text(strain_r(ind_stop(ii))+0, stress_r(ind_stop(ii))-10,['',num2str(ii)],'fontsize',18);
end
plot(0,0,'.r','markersize',18);
text(0,-10,['0'],'fontsize',18);

set(gca,'xdir','normal','ydir','normal','xlim',[-0.03 0.005],'ylim',[-150, 150]);
xlabel('Strain');ylabel('Stress (MPa)');
set(gca,'fontsize',18);

% legned1 = legend('Loading Interrupted Strain Level','location','east');
% set(legned1, 'Position',[0.24 0.30 0.65 0.075]);
%%
print(fullfile(output_dir, 'stress vs gage strain.tif'),'-dtiff');

%% save data
save(fullfile(output_dir, 'Mg4Al_S1_compression_curve_data'),'A','elongation_r','ind_stop','Load_r','strain_r','stress_r','time_r','ue_r');


%% using DIC strain to interpolate
exx_dic = [-0.0012   -0.0117   -0.0186   -0.0172   -0.0124    0.0006];

% find critical elongation and strainGage reading at the seven stops
elongation_at_dic = [];
ue_at_dic = [];
for ii=1:6
    elongation_at_dic(ii) = elongation(inds{ii}(end));
    ue_at_dic(ii) = ue(inds{ii}(end));
end

exx_dic = [0, exx_dic];
elongation_at_dic = [0, elongation_at_dic];
ue_at_dic = [0, ue_at_dic];

% Need to interpolate for each segment
exx_interped_by_dic = [];   % interp from displacement + DIC strain at load steps 
exx_interped_by_dic2 = [];  % interp from gage strain + DIC strain at load steps 
for ii = 1:6
    if ii == 1
        exx_interped_by_dic = [exx_interped_by_dic; interp1(elongation_at_dic(ii:ii+1), exx_dic(ii:ii+1), elongation(1:inds{ii}(end)),'linear','extrap')];
        exx_interped_by_dic2 = [exx_interped_by_dic2; interp1(ue_at_dic(ii:ii+1), exx_dic(ii:ii+1), ue(1:inds{ii}(end)),'linear','extrap')];
    else
        exx_interped_by_dic = [exx_interped_by_dic; interp1(elongation_at_dic(ii:ii+1), exx_dic(ii:ii+1), elongation(inds{ii-1}(end):inds{ii}(end)),'linear','extrap')];
        exx_interped_by_dic2 = [exx_interped_by_dic2; interp1(ue_at_dic(ii:ii+1), exx_dic(ii:ii+1), ue(inds{ii-1}(end):inds{ii}(end)),'linear','extrap')];
    end
end
% exx_interped_by_dic = interp1(elongation_at_dic, exx_dic, elongation);
exx_interped_by_dic_r = exx_interped_by_dic(rows,:);
exx_interped_by_dic_r2 = exx_interped_by_dic2(rows,:);

%% [plot 2] stress vs strain by displacement + DIC at load step
figure; hold on
colors = lines(6);

plot(exx_interped_by_dic_r, stress_r, '-', 'linewidth',2);
plot(0,0,'.r','markersize',18);
text(0, 0, ['0'],'fontsize',18);
for ii = 1:length(ind_stop)
   plot(exx_interped_by_dic_r(ind_stop(ii)),stress_r(ind_stop(ii)),'.r','markersize',18);
   text(exx_interped_by_dic_r(ind_stop(ii))+0, stress_r(ind_stop(ii))-10,['',num2str(ii)],'fontsize',18);
end

set(gca,'xdir','normal','ydir','normal','xlim',[-0.022 0.003],'ylim',[-150, 150]);
set(gca,'xTick', -0.02:0.005:0, 'XTickLabel',{'-0.02','-0.015','-0.01','-0.005','0'});

xlabel('Strain');ylabel('Stress (MPa)');
set(gca,'fontsize',18)

print(fullfile(output_dir, 'stress vs strain_by_disp_and_dic.tif'),'-dtiff');

%% [plot 3] stress vs strain by gage_strain + DIC at load step
figure; hold on
colors = lines(6);

plot(exx_interped_by_dic_r2, stress_r, '-', 'linewidth',2);
plot(0,0,'.r','markersize',18);
text(0, 0, ['0'],'fontsize',18);
for ii = 1:length(ind_stop)
   plot(exx_interped_by_dic_r2(ind_stop(ii)),stress_r(ind_stop(ii)),'.r','markersize',18);
   text(exx_interped_by_dic_r2(ind_stop(ii))+0, stress_r(ind_stop(ii))-10,['',num2str(ii)],'fontsize',18);
end

set(gca,'xdir','normal','ydir','normal','xlim',[-0.022 0.003],'ylim',[-150, 150]);
set(gca,'xTick', -0.02:0.005:0, 'XTickLabel',{'-0.02','-0.015','-0.01','-0.005','0'});

xlabel('Strain');ylabel('Stress (MPa)');
set(gca,'fontsize',18);

print(fullfile(output_dir, 'stress vs strain_by_gage_and_dic.tif'),'-dtiff');

%% display strain and stress values at intervals
t = table;
t.r1 = strain_r(ind_stop);
t.r2 = stress_r(ind_stop);
t.r3 = exx_interped_by_dic_r(ind_stop);
t.r4 = exx_interped_by_dic_r2(ind_stop);
t.Properties.VariableNames = {'sg_strain','stress','disp_dic_strain', 'gage_dic_strain'};
display(t);
%% plot up to strain-5
% figure; hold on
% 
% plot(0,0,'.r','markersize',18);
% text(-0.001,-18,['0'],'fontsize',18);
% 
% % plot up to strain level-5
% plot(exx_interped_by_dic_r(1:ind_stop(5)),stress_r(1:ind_stop(5)),'linewidth',2);
% 
% for ii = 1:length(ind_stop)-2
%    plot(exx_interped_by_dic_r(ind_stop(ii)),stress_r(ind_stop(ii)),'.r','markersize',18);
%    text(exx_interped_by_dic_r(ind_stop(ii))+0.001,stress_r(ind_stop(ii))-18,['',num2str(ii)],'fontsize',18);
% end
% 
% set(gca,'xdir','reverse','ydir','reverse','xlim',[-0.04 0],'ylim',[-250,0]);
% xlabel('Strain');ylabel('Stress (MPa)');
% set(gca,'fontsize',18);
%%
% print('stress strain curve 4.tif','-dtiff')


%%
save(fullfile(output_dir, 'Mg4Al_S1_compression_curve_data.mat'),'exx_interped_by_dic_r','exx_interped_by_dic_r2','exx_dic','-append');


%% get slop of elastic part
xx = strain_r(1:ind_stop(1));
yy = stress_r(1:ind_stop(1));
close all;
figure; hold on;
plot(0,0,'.r');
plot(xx,yy,'-');
[~,ind] = min(abs(yy-(-68)));
xxx = xx(1:ind);
yyy = yy(1:ind);
mdl = fitlm(xxx, yyy, 'y~x1-1');
elasticModulusExperimental = mdl.Coefficients.Estimate(1)/1000;
hold on;
drawline(gca, 'Position', [0, 0; xxx(end), xxx(end)*elasticModulusExperimental*1000], 'Color','r');
legend('Elastic Modulus by Strain Gage = 38.5 GPa','fontsize',12,'location','northwest');
print(fullfile(output_dir, 'elastic modulus by strain gage.tiff'),'-dtiff');


%% Look at what happend before load stop #1
close all;
% (1) How they change with time 
figure;
n = 173;

subplot(3,1,1);
plot(t_dds(1:n), elongation(1:n), '-or');
set(gca,'xlim',[0 190]);
ylabel('elongation, mm');
title('from 0 to load stop #1');

subplot(3,1,2);
plot(t_sg(1:n)-(t_sg(1)-t_dds(1)), ue(1:n), '-ob');
set(gca,'xlim',[0 190]);
ylabel('strain gage \mu\epsilon');

subplot(3,1,3);
plot(t_dds(1:n), stress(1:n), '-ok');
set(gca,'xlim',[0 190]);
ylabel('stress, MPa');
xlabel('Time, s')

% (2) How they change with stress
figure; 
subplot(1,3,1);
plot(elongation_r(1:n),stress(1:n),'-or');
set(gca,'xlim',[elongation_r(n), 0]);
xlabel('elongation,mm');
ylabel('stress, MPa');
title('from 0 to load stop #1');

subplot(1,3,2);
plot(exx_interped_by_dic_r(1:n),stress(1:n),'-or');
set(gca,'xlim',[exx_interped_by_dic_r(n), 0]);
xlabel('elongation to DIC strain');
ylabel('stress, MPa');

subplot(1,3,3);
plot(ue_r(1:n),stress(1:n),'-ob');
set(gca,'xlim',[ue_r(n), 0]);
xlabel('strain gage, \mu\epsilon');
ylabel('stress, MPa');

% can adjust and print








