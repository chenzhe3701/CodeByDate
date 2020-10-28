%
cd('C:\Users\ZheChen\Desktop\UMich work\Labview for Tensile Stage');

filename = 'Mg4Al_compression_2020_10_13.lvm';
delimiterIn = '\t';
headerlinesIn = 1;
A = importdata(filename, delimiterIn, headerlinesIn);
display(A.colheaders)

% [~, ind] = max(A.data(:,3));
% A.data = A.data(1:ind,:);
data = A.data;

%% make variable
time = A.data(:,1);
displacement =  A.data(:,2);
force =  A.data(:,3);
motor_speed = A.data(:,4);
motor_current = A.data(:,5);
strain =  A.data(:,6)/1000000;
stress = force/3.3/2.74;   % Mg4Al tested on 2020-10-13

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
figure;
plot(time, motor_speed);
clc;
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
ind_stop = [];
% load 0 -> 1
data = [data; A.data(5:1372,:)];
data = [data; A.data(3023:3172,:)];
ind_stop = [ind_stop, size(data,1)];
% load 1 -> 2
data = [data; A.data(7224:7399,:)];
ind_stop = [ind_stop, size(data,1)];
% load 2 -> 3
data = [data; A.data(11340:11518,:)];
ind_stop = [ind_stop, size(data,1)];
% load 3 -> 4
data = [data; A.data(15343:15687,:)];
ind_stop = [ind_stop, size(data,1)];

%% make variable, new
time = data(:,1);
displacement =  data(:,2);
force =  data(:,3);
motor_speed = data(:,4);
motor_current = data(:,5);
strain =  data(:,6)/1000000;
stress = force/3.3/2.74;   
%%
figure; hold on;
plot(strain, stress, 'linewidth', 3);
plot(strain(ind_stop(1:end-1)), stress(ind_stop(1:end-1)), 'r.', 'markersize', 18);
xlabel('Strain, from strain gage');
ylabel('Stress (MPa)');
set(gca, 'xlim',[-0.025, 0.001], 'ylim',[-120,0], 'fontsize',18);

figure; hold on;
plot(displacement, stress, 'linewidth', 3);
plot(displacement(ind_stop(1:end-1)), stress(ind_stop(1:end-1)), 'r.', 'markersize', 18);
xlabel('Displacement (mm)');
ylabel('Stress (MPa)');
set(gca, 'xlim',[-2, 0], 'ylim',[-120,0], 'fontsize',18);

