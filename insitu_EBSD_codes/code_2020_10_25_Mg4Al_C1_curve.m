% Analyze compression data for Mg4Al_C1, tested on 2020-10-23
close all;
clc;
sample_name = 'Mg4Al_C1';
working_dir = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu curve';
cd(working_dir);

filename = '2020-10-23 Mg4Al comp_ten_data.lvm';
delimiterIn = '\t';
headerlinesIn = 1;
A = importdata(filename, delimiterIn, headerlinesIn);
display(A.colheaders);
disp('size of A.data:')
disp(size(A.data,1));

filename_2 = '2020-10-23 Mg4Al comp_ten_final_to_0N.lvm';
B = importdata(filename_2, delimiterIn, headerlinesIn);
B.data(:,1) = B.data(:,1) + A.data(end,1) + 1;  % remake B time

% [~, ind] = max(A.data(:,3));
% A.data = A.data(1:ind,:);
A.data = [A.data; B. data]; % combine B.data into A.data
A.data(:,6) = A.data(:,6) - A.data(1,6);    % zero strain

width = 3.48;       % Mg4Al_C1 tested on 2020-10-13
thickness = 2.69;
%% make variable
time = A.data(:,1);
displacement =  A.data(:,2);
force =  A.data(:,3);
motor_speed = A.data(:,4);
motor_current = A.data(:,5);
strain =  A.data(:,6)/1000000;
stress = force/width/thickness;   

%% explore
close all;

figure;
plot(time, displacement, '-r', 'linewidth', 3);
xlabel('time (sec)');
ylabel('displacement (mm)');
set(gca,'fontsize',16);

figure;
plot(time, force, '-r', 'linewidth', 3);
xlabel('time (sec)');
ylabel('load (N)');
set(gca,'fontsize',16);

figure;
plot(time, strain, '-r', 'linewidth', 3);
xlabel('time (sec)');
ylabel('strain');
set(gca,'fontsize',16);

%% Select active loading part, exploratory
disp('Select loading parts:');
figure;
plot(time, motor_speed);
indn_max = length(motor_speed);
ind = 1;
while ind < indn_max
    % find starting point of a segment with motor speed > 0.002
    if abs(motor_speed(ind)) < 0.003
        ind = ind + 1;
    else
        ind_a = ind;    % starting point

        % find ending point of this segment
        while (ind < indn_max) && (abs(motor_speed(ind)) > 0.003)
            ind_b = ind;    % ending point
            ind = ind + 1;
        end
        disp([ind_a, ind_b]);
    end
end
%% Construct active loading part
data = [];
ind_stop = 1;
% load 0 -> 1
data = [data; A.data(5:781,:)];
ind_stop = [ind_stop, size(data,1)];
% load 1 -> 2
data = [data; A.data(10809:10958,:)];
ind_stop = [ind_stop, size(data,1)];
% load 2 -> 3
data = [data; A.data(21122:21309,:)];
ind_stop = [ind_stop, size(data,1)];
% load 3 -> 4
data = [data; A.data(32195:32567,:)];
ind_stop = [ind_stop, size(data,1)];

% load 4 -> 5
data = [data; A.data(42738:44288,:)];
ind_stop = [ind_stop, size(data,1)];
% load 5 -> 6
data = [data; A.data(54502:55452,:)];
ind_stop = [ind_stop, size(data,1)];
% load 6 -> 7
data = [data; A.data(65108:65692,:)];
ind_stop = [ind_stop, size(data,1)];
% load 7 -> 8
data = [data; A.data(65726:66207,:)];
ind_stop = [ind_stop, size(data,1)];

%% make variable, new, and subplot [disp,force,strain] vs. [time]
time = data(:,1);
displacement =  data(:,2);
force =  data(:,3);
motor_speed = data(:,4);
motor_current = data(:,5);
strain =  data(:,6)/1000000;
stress = force/width/thickness;   

figure; set(gcf,'Position', [150, 150, 600, 800]);
subplot(3,1,1); hold on;
plot(displacement, '-r', 'linewidth', 3);
plot(ind_stop(1:end-1), displacement(ind_stop(1:end-1)), '.b', 'markersize',16);
xlabel('time (sec)');
ylabel('displacement (mm)');
set(gca,'fontsize',12);

subplot(3,1,2); hold on;
plot(force, '-r', 'linewidth', 3);
plot(ind_stop(1:end-1), force(ind_stop(1:end-1)), '.b', 'markersize',16);
xlabel('time (sec)');
ylabel('load (N)');
set(gca,'fontsize',12);

subplot(3,1,3); hold on;
plot(strain, '-r', 'linewidth', 3);
plot(ind_stop(1:end-1), strain(ind_stop(1:end-1)), '.b', 'markersize',16);
xlabel('time (sec)');
ylabel('strain');
set(gca,'fontsize',12);

disp('stress at load steps:');
disp(stress(ind_stop));
disp('strain at load steps:');
disp(strain(ind_stop));
disp('displacement at load steps:');
disp(displacement(ind_stop));
print(fullfile(working_dir,'displacement load strain vs time.tiff'),'-dtiff');
%% stress vs strain,  displacement vs strain
figure; hold on;
% plot(strain, stress, 'linewidth', 3);
colors = parula(length(ind_stop));
for ii = 1:length(ind_stop)-1
    plot(strain(ind_stop(ii):ind_stop(ii+1)), stress(ind_stop(ii):ind_stop(ii+1)), 'color',colors(ii,:), 'linewidth',3);
end
plot(strain(ind_stop(1:end-1)), stress(ind_stop(1:end-1)), 'r.', 'markersize', 18);
xlabel('Strain, from strain gage');
ylabel('Stress (MPa)');
set(gca, 'xlim',[-0.028, 0.001], 'ylim',[-150,150], 'fontsize',18);
print(fullfile(working_dir,'stress vs strain.tiff'),'-dtiff');

figure; hold on;
% plot(displacement, stress, 'linewidth', 3);
for ii = 1:length(ind_stop)-1
    plot(displacement(ind_stop(ii):ind_stop(ii+1)), stress(ind_stop(ii):ind_stop(ii+1)), 'color',colors(ii,:), 'linewidth',3);
end
plot(displacement(ind_stop(1:end-1)), stress(ind_stop(1:end-1)), 'r.', 'markersize', 18);
xlabel('Displacement (mm)');
ylabel('Stress (MPa)');
set(gca, 'xlim',[-1.5, 1.5], 'ylim',[-150,150], 'fontsize',18);
print(fullfile(working_dir,'stress vs displacement.tiff'),'-dtiff');

tbl_full = array2table([stress(:),strain(:),displacement(:)],'VariableNames',{'stress','strain_sg','displacement'});
tbl = array2table([[0:length(ind_stop)-1]',stress(ind_stop),strain(ind_stop),displacement(ind_stop)],'VariableNames',{'iE','stress','strain_sg','displacement'});
disp(tbl);
figure;
uitable('Data',tbl{:,:},'ColumnName',tbl.Properties.VariableNames,...
    'RowName',tbl.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
print(fullfile(working_dir,'stress strain table.tiff'),'-dtiff');

save(fullfile(working_dir, [sample_name,'_processed_loading_data.mat']), 'displacement','stress','strain','ind_stop','tbl');

%% [] stress vs strain, each load step with circle indicating position
clc;
for iE = 0:8
    figure; hold on;
    colors = parula(15);
    for ii = 1:8
        plot(strain(ind_stop(ii):ind_stop(ii+1)), stress(ind_stop(ii):ind_stop(ii+1)), 'color',colors(ii,:), 'linewidth',3);
    end
    plot(strain(ind_stop(1:9)), stress(ind_stop(1:9)), 'r.', 'markersize', 18);
    plot(strain(ind_stop(iE+1)), stress(ind_stop(iE+1)), 'ro', 'markersize', 12, 'linewidth',3);
    xlabel('Strain, from strain gage');
    ylabel('Stress (MPa)');
    set(gca, 'xlim',[-0.028, 0.001], 'ylim',[-150,150], 'fontsize',18);
    print(fullfile(working_dir,['stress vs strain iE_',num2str(iE),'.tiff']),'-dtiff');
    close;
end
%%
close all;

