
%% (1) test drift of stage. At controlled displacement rate of 0, the displacement, load, and strain reading kept changing overnight.

close all;
clc;

working_dir = 'E:\zhec umich Drive\2020-11-13 Mg4Al_C2 insiu curve';
cd(working_dir);

fileName = '2020-11-13 test_drift.lvm';
delimiterIn = '\t';
headerlinesIn = 1;
A = importdata(fullfile(working_dir,fileName), delimiterIn, headerlinesIn);
display(A);

%% make variable
time = A.data(:,1);
displacement =  A.data(:,2);
force =  A.data(:,3);
motor_speed = A.data(:,4);
motor_current = A.data(:,5);
strain =  A.data(:,6)/1000000;
stress = force/3.44/2.68;   % Mg4Al_C2 tested on 2020-11-13

%% explore drift during drift test
close all;

figure;
subplot(3,1,1); 
plot(time, displacement, '-r', 'linewidth', 3);
xlabel('time (sec)');
ylabel('displacement (mm)');
set(gca,'fontsize',12);

subplot(3,1,2); 
plot(time, force, '-r', 'linewidth', 3);
xlabel('time (sec)');
ylabel('load (N)');
set(gca,'fontsize',12);

subplot(3,1,3); 
plot(time, strain, '-r', 'linewidth', 3);
xlabel('time (sec)');
ylabel('strain');
set(gca,'fontsize',12);

% figure;
% plot(time, displacement);
% set(gca,'xlim',[0 1800], 'ylim', [0.395, 0.405]);
mdl = fitlm(time, displacement);
drift_rate_micron_per_hour = mdl.Coefficients.Estimate(2) * 1000 * 3600
mdl = fitlm(time, force);
drift_rate_newton_per_hour = mdl.Coefficients.Estimate(2) * 3600
mdl = fitlm(time, strain);
drift_rate_strain_per_hour = mdl.Coefficients.Estimate(2) * 3600

%% (2) Look at the drift during a real test

fileName = '2020-11-13 Mg4Al_C2 comp_ten_data.lvm';

A = importdata(fullfile(working_dir,fileName), delimiterIn, headerlinesIn);
display(A);

%% make variable
time = A.data(:,1);
displacement =  A.data(:,2);
force =  A.data(:,3);
motor_speed = A.data(:,4);
motor_current = A.data(:,5);
strain =  A.data(:,6)/1000000;
stress = force/3.44/2.68;   % Mg4Al_C2 tested on 2020-11-13

%% explore
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

%% Show active loading part, exploratory
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

%% Select the part where displacement was hold
close; close; close; close
data = A.data(1484:end,:);

%% make variable, new, and subplot [disp,force,strain] vs. [time]
time = data(:,1) - data(1,1);
displacement =  data(:,2);
force =  data(:,3);
motor_speed = data(:,4);
motor_current = data(:,5);
strain =  data(:,6)/1000000;
stress = force/3.3/2.74;   

figure;
subplot(3,1,1); hold on;
plot(time, displacement, '-r', 'linewidth', 3);
xlabel('time (sec)');
ylabel('displacement (mm)');
set(gca,'fontsize',12);

subplot(3,1,2); hold on;
plot(time, force, '-r', 'linewidth', 3);
xlabel('time (sec)');
ylabel('load (N)');
set(gca,'fontsize',12);

subplot(3,1,3); hold on;
plot(time, strain, '-r', 'linewidth', 3);
xlabel('time (sec)');
ylabel('strain');
set(gca,'fontsize',12);

% figure;
% plot(time, displacement);
% set(gca,'xlim',[0 1800], 'ylim', [-1.395, -1.385]);
mdl = fitlm(time, displacement);
drift_rate_micron_per_hour = mdl.Coefficients.Estimate(2) * 1000 * 3600
mdl = fitlm(time, force);
drift_rate_newton_per_hour = mdl.Coefficients.Estimate(2) * 3600
mdl = fitlm(time, strain);
drift_rate_strain_per_hour = mdl.Coefficients.Estimate(2) * 3600