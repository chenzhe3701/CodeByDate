% Look at the stress-strain curve of WE43-T6-C1 again.
% 2021-01-10

%% strain gage data
clear;clc;
close all;
data_dir = 'D:\WE43_T6_C1\Compression Data';
save_dir = 'C:\Users\ZheChen\Desktop';

fid = fopen(fullfile(data_dir, 'WE43_T6_C1_compression_strain_gage.TXT'));
a = textscan(fid,'%{hh:mm:ss}D %f','delimiter','\t','headerlines',6); % determined that it started 9sec earlier
ue = a{1,2};
t_sg = (0:length(ue)-1)';
fclose(fid);

ind_reset_sg = find(ue==-34119);     % reset at this point, when ue = -34119
ue(ind_reset_sg+1:end) = ue(ind_reset_sg+1:end) + ue(ind_reset_sg);

%% DDS data
fid = fopen(fullfile(data_dir, 'WE43_T6_C1_insitu_compression.TXT'));
b = textscan(fid,'%s %f %f %f %f %f','delimiter','\t','headerlines',2);
Load = b{1,3};
elongation = b{1,4};
t_dds = (0:length(Load)-1)';
fclose(fid);

A = 10.9;
stress = Load/A;

%% use for plot raw data
figure;
plot(t_sg,ue/min(ue(:)),'r');
hold;
plot(t_dds,elongation/min(elongation(:)),'b');
legend({'ue','elongation'});

%% make them the same length (device time runs in different speed ...). Note: both aligned to [dds time].
ind_ue = (8:187986)';
ind_dds = (1:187992)';

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
%% select data for active loading part
inds = [{13:363}, {29805:29854},  {[56796:56827,57292:57306,58144:58168]},...
    {84932:85074}, {112071:112243}, {136710:137081}, {161593:162218}];
rs = [];
for ii = 1:length(inds)
    rs = [rs,inds{ii}];
    ind_stop(ii) = length(rs);
end

% rs = [13:363, 29805:29854,  56796:56827,57292:57306,58144:58168,...
%     84932:85074, 112071:112243, 136710:137081, 161593:162218];

% reduced
time_r = t_dds(rs,:);
Load_r = Load(rs,:);
stress_r = stress(rs,:);
elongation_r = elongation(rs,:);
ue_r = ue(rs,:);

figure; plot(elongation_r, ue_r);

[temp,idx] = min(abs(ue_r-(-34119)));

%% interp/extrap ue with elongation
inds_interp = 260:690;
model = fitlm(elongation_r(inds_interp,:),ue_r(inds_interp,:));
slope = model.Coefficients.Estimate(2);
ue_r(inds_interp(end)+1:end) = elongation_r(inds_interp(end)+1:end) * slope;
strain_r = ue_r/1000000;
%% plot
figure; hold on

plot(0,0,'.r','markersize',18);
text(-0.002,-18,['0'],'fontsize',18);
plot(strain_r,stress_r,'linewidth',2);
for ii = 1:length(ind_stop)
   plot(strain_r(ind_stop(ii)),stress_r(ind_stop(ii)),'.r','markersize',18);
   text(strain_r(ind_stop(ii))+0.003,stress_r(ind_stop(ii))-18,['',num2str(ii)],'fontsize',18);
end

set(gca,'xdir','reverse','ydir','reverse','xlim',[-0.065 0]);
xlabel('Strain');ylabel('Stress (MPa)');
set(gca,'fontsize',18)

% legned1 = legend('Loading Interrupted Strain Level','location','east');
% set(legned1, 'Position',[0.24 0.30 0.65 0.075]);
%%
print('stress strain curve 2.tif','-dtiff')

%% save data
save('compression_curve_data','A','elongation_r','ind_stop','Load_r','strain_r','stress_r','time_r','ue_r');


%% using DIC strain to interpolate
exx_dic = [-0.003, -0.0036, -0.0115, -0.023, -0.0385, -0.0567, -0.0829];

% find critical elongation and strainGage reading at the seven stops
elongation_at_dic = [];
ue_at_dic = [];
for ii=1:7
    elongation_at_dic(ii) = elongation(inds{ii}(end));
    ue_at_dic(ii) = ue(inds{ii}(end));
end

exx_dic = [0, exx_dic];
elongation_at_dic = [0, elongation_at_dic];

exx_interped_by_dic = interp1(elongation_at_dic, exx_dic, elongation);

exx_interped_by_dic_r = exx_interped_by_dic(rs,:);

%% plot
figure; hold on

plot(0,0,'.r','markersize',18);
text(-0.002,-18,['0'],'fontsize',18);
plot(exx_interped_by_dic_r,stress_r,'linewidth',2);

for ii = 1:length(ind_stop)
   plot(exx_interped_by_dic_r(ind_stop(ii)),stress_r(ind_stop(ii)),'.r','markersize',18);
   text(exx_interped_by_dic_r(ind_stop(ii))+0.003,stress_r(ind_stop(ii))-18,['',num2str(ii)],'fontsize',18);
end

set(gca,'xdir','reverse','ydir','reverse','xlim',[-0.09 0],'ylim',[-300,0]);
xlabel('Strain');ylabel('Stress (MPa)');
set(gca,'fontsize',18);
%%
print('stress strain curve 3.tif','-dtiff')

% legned1 = legend('Loading Interrupted Strain Level','location','east');
% set(legned1, 'Position',[0.24 0.30 0.65 0.075]);

%% plot up to strain-5, for 2021 Twin Stats paper
figure; hold on

plot(0,0,'.r','markersize',18);
text(-0.001,-18,['0'],'fontsize',18);

% plot up to strain level-5
plot(exx_interped_by_dic_r(1:ind_stop(5)),stress_r(1:ind_stop(5)),'linewidth',2,'color','b');

for ii = 1:length(ind_stop)-2
   plot(exx_interped_by_dic_r(ind_stop(ii)),stress_r(ind_stop(ii)),'.r','markersize',18);
   text(exx_interped_by_dic_r(ind_stop(ii))+0.001,stress_r(ind_stop(ii))-18,['',num2str(ii)],'fontsize',18);
end

set(gca,'xdir','reverse','ydir','reverse','xlim',[-0.04 0],'ylim',[-250,0]);
xlabel('Strain');ylabel('Stress (MPa)');
set(gca,'fontsize',18);
%%
print(fullfile(save_dir, 'WE43_T6 curve TwinStats paper 2021.tif'),'-dtiff');


%%
save('compression_curve_data','exx_interped_by_dic_r','exx_dic','-append');


%% Plot strain levels 0,2,3,4,5 for 2020 Mater Char paper
%% plot up to strain-5
figure; hold on

plot(0,0,'.r','markersize',18);
text(-0.001,-18,['1'],'fontsize',18);

% plot up to strain level-5
plot(exx_interped_by_dic_r(1:ind_stop(5)),stress_r(1:ind_stop(5)),'linewidth',2,'color','b');

for ii = 2:5
   plot(exx_interped_by_dic_r(ind_stop(ii)),stress_r(ind_stop(ii)),'.r','markersize',18);
   text(exx_interped_by_dic_r(ind_stop(ii))+0.001,stress_r(ind_stop(ii))-18,['',num2str(ii)],'fontsize',18);
end

set(gca,'xdir','reverse','ydir','reverse','xlim',[-0.04 0],'ylim',[-250,0]);
xlabel('Strain');ylabel('Stress (MPa)');
set(gca,'fontsize',18);

%%
print(fullfile(save_dir, 'WE43_T6 curve for MaterChar 2020.tif'),'-r300','-dtiff');



