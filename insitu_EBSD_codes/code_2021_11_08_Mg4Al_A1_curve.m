%% load data

close all;
clc;
sample_name = 'Mg4Al_A1';
working_dir = 'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu curve';
cd(working_dir);

fileName = '2021-10-28 Mg4Al_A1 compression-tension.lvm';
delimiterIn = '\t';
headerlinesIn = 1;
A = importdata(fullfile(working_dir,fileName), delimiterIn, headerlinesIn);
display(A);

width = 2.74;       % Mg4Al_A1, tested on 2021-10-28
thickness = 2.86;

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
speed_sign = [];
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
        % positive motor_speed signal = compression
        speed_sign = [speed_sign; sign(mean(motor_speed(ind_a:ind_b)))];
        disp([num2str(ind_a),',',num2str(ind_b),';']);
    end
end
disp(speed_sign);


%% Construct active loading part 

inds = [5,591;
617,625;    % 0->1
8206,8397;
8414,8424;  % 1->2
15779,16024;
16042,16052;    % 2->3
23562,23687;    % 3->4
31069,32821;
32836,32843;    % 4->5
40243,40559;
40578,40591;    % 5->6
49017,49274;
49290,49326;    % 6->7
57389,59464;
59484,59493;    % 7->8
66952,67157;
67173,67183;    % 8->9
74554,74683;
74685,74686;    
74689,74690;
74694,74694;    
74699,74699;
74704,74852;    
74866,74881;    % 9->10
82567,83791;
83809,83818;    % 10->11
91265,91568;
91583,91608;    % 11->12
99004,99262;
99277,99316;    % 12->13
107056,107309;  % 13->14
];

% is the index a stop, where in-situ test was paused for imaging?
is_stop = [0,1, 0,1, 0,1, 1, 0,1,   0,1, 0,1, 0,1, 0,1, 0,0,0,0,0,0,1,  0,1, 0,1, 0,1, 1]; 

data = [];
ind_stop = 1;

for ir = 1:size(inds,1)
    data = [data; A.data(inds(ir,1):inds(ir,2),:)];
    if is_stop(ir)
        ind_stop = [ind_stop, size(data,1)];
    end
    % max iE that we want to include
    max_iE_to_include = 16;
    if length(ind_stop)>max_iE_to_include
        break;
    end
end


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
set(gca,'fontsize',12,'xlim',[0 10000]);

subplot(3,1,2); hold on;
plot(force, '-r', 'linewidth', 3);
plot(ind_stop(1:end-1), force(ind_stop(1:end-1)), '.b', 'markersize',16);
xlabel('time (sec)');
ylabel('load (N)');
set(gca,'fontsize',12,'xlim',[0 10000]);

subplot(3,1,3); hold on;
plot(strain, '-r', 'linewidth', 3);
plot(ind_stop(1:end-1), strain(ind_stop(1:end-1)), '.b', 'markersize',16);
xlabel('time (sec)');
ylabel('strain');
set(gca,'fontsize',12,'xlim',[0 10000]);

disp('stress at load steps:');
disp(stress(ind_stop));
disp('strain at load steps:');
disp(strain(ind_stop));
disp('displacement at load steps:');
disp(displacement(ind_stop));
print(fullfile(working_dir,'displacement load strain vs time.tiff'),'-dtiff');
%% stress vs strain,  displacement vs strain
% [] stress vs strain gage strain
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
print(fullfile(working_dir,'stress vs strain.tiff'),'-dtiff');

% [] stress vs displacement
figure; hold on;
for ii = 1:length(ind_stop)-1
    plot(displacement(ind_stop(ii):ind_stop(ii+1)), stress(ind_stop(ii):ind_stop(ii+1)), 'color',colors(ii,:), 'linewidth',3);
end
plot(displacement(ind_stop(1:end-1)), stress(ind_stop(1:end-1)), 'r.', 'markersize', 18);
xlabel('Displacement (mm)');
ylabel('Stress (MPa)');
set(gca, 'xlim',[-1.3, 1.6], 'ylim',[-150,150], 'fontsize',18);
print(fullfile(working_dir,'stress vs displacement.tiff'),'-dtiff');

% [] table summary
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
for iE = 0:13
    figure; hold on;
    colors = parula(15);
    for ii = 1:13
        plot(strain(ind_stop(ii):ind_stop(ii+1)), stress(ind_stop(ii):ind_stop(ii+1)), 'color',colors(ii,:), 'linewidth',3);
    end
    plot(strain(ind_stop(1:14)), stress(ind_stop(1:14)), 'r.', 'markersize', 18);
    plot(strain(ind_stop(iE+1)), stress(ind_stop(iE+1)), 'ro', 'markersize', 12, 'linewidth',3);
    xlabel('Strain, from strain gage');
    ylabel('Stress (MPa)');
    set(gca, 'xlim',[-0.03, 0.005], 'ylim',[-150,150], 'fontsize',18);
    print(fullfile(working_dir,['stress vs strain iE_',num2str(iE),'.tiff']),'-dtiff');
    close;
end
%%
close all;



