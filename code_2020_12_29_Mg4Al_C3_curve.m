% Analyze compression data for UM_134_Mg_C1, tested on 2020-12-05
%% load data

close all;
clc;

working_dir = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu curve';
cd(working_dir);

fileName = '2020-12-23 Mg4Al_C3 comp_ten_data.lvm';
delimiterIn = '\t';
headerlinesIn = 1;
A = importdata(fullfile(working_dir,fileName), delimiterIn, headerlinesIn);
display(A);

width = 3.47;  % Mg4Al_C3 tested on 2020-12-23
thickness = 2.64;
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
data = [data; A.data(4:740,:)];
ind_stop = [ind_stop, size(data,1)];
% load 1 -> 2
data = [data; A.data(10795:10992,:)];
ind_stop = [ind_stop, size(data,1)];
% load 2 -> 3
data = [data; A.data(21206:21478,:)];
ind_stop = [ind_stop, size(data,1)];
% load 3 -> 4
data = [data; A.data(31270:31421,:)];
ind_stop = [ind_stop, size(data,1)];
% load 4 -> 5
data = [data; A.data(41242:43007,:)];
ind_stop = [ind_stop, size(data,1)];
% load 5 -> 6
data = [data; A.data(52798:53159,:)];
ind_stop = [ind_stop, size(data,1)];
% load 6 -> 7
data = [data; A.data(68899:69219,:)];
ind_stop = [ind_stop, size(data,1)];
% load 7 -> 8
data = [data; A.data(79018:81262,:)];
ind_stop = [ind_stop, size(data,1)];
% load 8 -> 9
data = [data; A.data(91744:91966,:)];
ind_stop = [ind_stop, size(data,1)];
% load 9 -> 10
data = [data; A.data(101777:102062,:)];
ind_stop = [ind_stop, size(data,1)];
% load 10 -> 11
data = [data; A.data(113578:115844,:)];
ind_stop = [ind_stop, size(data,1)];
% load 11 -> 12
data = [data; A.data(125821:126133,:)];
ind_stop = [ind_stop, size(data,1)];
% load 12 -> 13
data = [data; A.data(135995:136257,:)];
ind_stop = [ind_stop, size(data,1)];
% load 13 -> 14
data = [data; A.data(146336:146694,:)];
ind_stop = [ind_stop, size(data,1)];
%% make variable, new, and subplot [disp,force,strain] vs. [time]
time = data(:,1);
displacement =  data(:,2);
force =  data(:,3);
motor_speed = data(:,4);
motor_current = data(:,5);
strain =  data(:,6)/1000000;
stress = force/width/thickness;

figure; set(gcf,'Position', [150, 150, 1024, 768]);
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
% plot(strain, stress, 'linewidth', 3);
colors = parula(length(ind_stop));
for ii = 1:length(ind_stop)-1
    plot(strain(ind_stop(ii):ind_stop(ii+1)), stress(ind_stop(ii):ind_stop(ii+1)), 'color',colors(ii,:), 'linewidth',3);
end
plot(strain(ind_stop(1:end-1)), stress(ind_stop(1:end-1)), 'r.', 'markersize', 18);
xlabel('Strain, from strain gage');
ylabel('Stress (MPa)');
set(gca, 'xlim',[-0.03, 0.005], 'ylim',[-150,150], 'fontsize',18);

figure; hold on;
for ii = 1:length(ind_stop)-1
    plot(displacement(ind_stop(ii):ind_stop(ii+1)), stress(ind_stop(ii):ind_stop(ii+1)), 'color',colors(ii,:), 'linewidth',3);
end
plot(displacement(ind_stop(1:end-1)), stress(ind_stop(1:end-1)), 'r.', 'markersize', 18);
xlabel('Displacement (mm)');
ylabel('Stress (MPa)');
set(gca, 'xlim',[-1.5, 1.5], 'ylim',[-150,150], 'fontsize',18);


%% Then, the strain was estimated from EBSD map, and corrected.

% from python, EBSD estimated strain

scale = [1, 0.9892, 0.9782, 0.9700, ...
    0.9704, 0.9759, 0.9880, 0.9975, ...
    0.9900, 0.9808, 0.9688, ...
    0.9775, 0.9876, 0.9964] ;

strain_ebsd = scale - 1;
strain_sg = strain(ind_stop(1:end-1));

figure; hold on;
for iE = 1:13
    plot(strain_sg(iE:iE+1), strain_ebsd(iE:iE+1), '-', 'color', colors(iE,:), 'linewidth', 1);
end
plot(strain_sg, strain_ebsd,'.r','markersize',24);

set(gca,'xlim',[-0.035,0.005], 'ylim',[-0.035,0.005], 'fontsize',18);
xlabel('Strain gage strain');
ylabel('EBSD estimated strain');












