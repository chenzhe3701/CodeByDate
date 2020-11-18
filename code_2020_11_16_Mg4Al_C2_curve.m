%% load data

close all;
clc;

working_dir = 'E:\zhec umich Drive\2020-11-13 Mg4Al_C2 insiu curve';
cd(working_dir);

fileName_1 = '2020-11-13 Mg4Al_C2 comp_ten_data.lvm';
delimiterIn = '\t';
headerlinesIn = 1;
A = importdata(fullfile(working_dir,fileName_1), delimiterIn, headerlinesIn);
display(A);

fileName_2 = '2020-11-13 Mg4Al_C2 comp_ten_data part_2.lvm';
delimiterIn = '\t';
headerlinesIn = 1;
B = importdata(fullfile(working_dir,fileName_1), delimiterIn, headerlinesIn);
display(B);

% concatenate B with A, for time (col_1), displacement (col_2), load (col_3), strain (col_6)
B.data(:,1) = B.data(:,1) - B.data(1,1) + A.data(end,1) + 1;
B.data(:,2) = B.data(:,2) - B.data(1,2) + A.data(end,2); 
B.data(:,3) = B.data(:,3) - B.data(1,3) + A.data(end,3); 
B.data(:,6) = B.data(:,6) - B.data(1,6) + A.data(end,6); 

A.data = [A.data;B.data];
%% make variable
time = A.data(:,1);
displacement =  A.data(:,2);
force =  A.data(:,3);
motor_speed = A.data(:,4);
motor_current = A.data(:,5);
strain =  A.data(:,6)/1000000;
stress = force/3.44/2.68;   % Mg4Al_C2 tested on 2020-11-13

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
data = [data; A.data(5:1483,:)];
ind_stop = [ind_stop, size(data,1)];
% load 1 -> 2
data = [data; A.data(3031:4509,:)];
ind_stop = [ind_stop, size(data,1)];

%% make variable, new, and subplot [disp,force,strain] vs. [time]
time = data(:,1);
displacement =  data(:,2);
force =  data(:,3);
motor_speed = data(:,4);
motor_current = data(:,5);
strain =  data(:,6)/1000000;
stress = force/3.3/2.74;   

figure;
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

%% stress vs strain,  displacement vs strain
figure; hold on;
plot(strain, stress, 'linewidth', 3);
plot(strain(ind_stop(1:end-1)), stress(ind_stop(1:end-1)), 'r.', 'markersize', 18);
xlabel('Strain, from strain gage');
ylabel('Stress (MPa)');
set(gca, 'xlim',[-0.02, 0.001], 'ylim',[-200,50], 'fontsize',18);

figure; hold on;
plot(displacement, stress, 'linewidth', 3);
plot(displacement(ind_stop(1:end-1)), stress(ind_stop(1:end-1)), 'r.', 'markersize', 18);
xlabel('Displacement (mm)');
ylabel('Stress (MPa)');
set(gca, 'xlim',[-3, 0.5], 'ylim',[-200,50], 'fontsize',18);


















