% look at weird relaxation data
% First, similar to analyzing loading curve, get regular active loading data  
working_dir = 'E:\zhec umich Drive\2021-02-03 test relaxation using UM134 Mg_C3';

fileName = '2020-02-03 test relaxation.lvm';
delimiterIn = '\t';
headerlinesIn = 1;
A = importdata(fullfile(working_dir,fileName), delimiterIn, headerlinesIn);
display(A);

width = 3.45;       % UM134_Mg_C3, tested on 2021-01-29
thickness = 2.63;
%% make variable
time = A.data(:,1);
displacement =  A.data(:,2);
force =  A.data(:,3);
motor_speed = A.data(:,4);
motor_current = A.data(:,5);
strain =  A.data(:,6)/1000000;
stress = force/width/thickness;   

%% explore. Note that the displacement sensor did not work, we need to manually fill displacement data ... 
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
    % find starting point of a segment with motor speed > 0.005
    if abs(motor_speed(ind)) < 0.005
        ind = ind + 1;
    else
        ind_a = ind;    % starting point

        % find ending point of this segment
        while (ind < indn_max) && (abs(motor_speed(ind)) > 0.005)
            ind_b = ind;    % ending point
            ind = ind + 1;
        end
        disp([num2str(ind_a),',',num2str(ind_b),';']);
    end
end

%% For this experiment, the displacement reading did not change well during test, so we need to manually fix
% for active loading part, the rate is 1 um/s
% for other parts, the rate is 0 um/s
inds = [4,784;
13268,13670];

sign_v = [1,-1]; % increase or decrease displacement  
disp_v = zeros(size(A.data,1),1);
for ir = 1:size(inds,1)
   disp_local = (1:(inds(ir,2)-inds(ir,1)+1)) * sign_v(ir);   % calculate local series
   disp_v(inds(ir,1):inds(ir,2)) = disp_v(inds(ir,1)) + disp_local; % add local increase
   disp_v(inds(ir,2):end) = disp_v(inds(ir,2));     % from pos_index_2, copy the value at pos_index_2
end

A.data(:,2) = disp_v / 1000;

%% Construct active loading part
data = [];
ind_stop = 1;
% load 0 -> 1
data = [data; A.data(4:784,:)];
ind_stop = [ind_stop, size(data,1)];
% load 1 -> 2
data = [data; A.data(13268:13670,:)];
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
set(gca,'fontsize',12,'xlim',[0 1300]);

subplot(3,1,2); hold on;
plot(force, '-r', 'linewidth', 3);
plot(ind_stop(1:end-1), force(ind_stop(1:end-1)), '.b', 'markersize',16);
xlabel('time (sec)');
ylabel('load (N)');
set(gca,'fontsize',12,'xlim',[0 1300]);

subplot(3,1,3); hold on;
plot(strain, '-r', 'linewidth', 3);
plot(ind_stop(1:end-1), strain(ind_stop(1:end-1)), '.b', 'markersize',16);
xlabel('time (sec)');
ylabel('strain');
set(gca,'fontsize',12,'xlim',[0 1300]);

disp('stress at load steps:');
disp(stress(ind_stop));
disp('strain at load steps:');
disp(strain(ind_stop));
disp('displacement at load steps:');
disp(displacement(ind_stop));

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
set(gca, 'xlim',[-0.005, 0.01], 'ylim',[-10,100], 'fontsize',18);
print(fullfile(working_dir,'stress vs strain.tiff'),'-dtiff');

figure; hold on;
for ii = 1:length(ind_stop)-1
    plot(displacement(ind_stop(ii):ind_stop(ii+1)), stress(ind_stop(ii):ind_stop(ii+1)), 'color',colors(ii,:), 'linewidth',3);
end
plot(displacement(ind_stop(1:end-1)), stress(ind_stop(1:end-1)), 'r.', 'markersize', 18);
xlabel('Displacement (mm)');
ylabel('Stress (MPa)');
set(gca, 'xlim',[-0.2, 1], 'ylim',[-10,100], 'fontsize',18);
print(fullfile(working_dir,'stress vs displacement.tiff'),'-dtiff');

