% look at the non-active/holding part of curve of UM134 Mg_C3

close all;
clc;
save_dir = 'E:\zhec umich Drive\0_temp_output\2021-03-22 relaxation analysis';
mkdir(save_dir);

%% [1] UM134 Mg_C3

working_dir = 'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu curve';
fileName = '2021-01-29 UM134_Mg_C3 comp_ten_data.lvm';
delimiterIn = '\t';
headerlinesIn = 1;
A = importdata(fullfile(working_dir,fileName), delimiterIn, headerlinesIn);
display(A);

width = 3.45;       % UM134_Mg_C3, tested on 2021-01-29
thickness = 2.63;

% For this experiment, the displacement reading did not change well during test, so we need to manually fix
% For active loading part, the rate is 1 um/s. For other parts, the rate is 0 um/    
inds = [4,1492;
11103,11226;
20903,21076;
31003,32266;
41903,42464;
52603,52809;
63103,63294;
73504,75454;
88804,88956;
114854,115036;
126003,127947;
138203,138397;
148703,148885;
159273,159407];

sign_v = [-1,-1,-1,1,1,1,1,-1,-1,-1,1,1,1,-1]; % increase or decrease displacement  
disp_v = zeros(size(A.data,1),1);
for ir = 1:size(inds,1)
   disp_local = (1:(inds(ir,2)-inds(ir,1)+1)) * sign_v(ir);   % calculate local series
   disp_v(inds(ir,1):inds(ir,2)) = disp_v(inds(ir,1)) + disp_local; % add local increase
   disp_v(inds(ir,2):end) = disp_v(inds(ir,2));     % from pos_index_2, copy the value at pos_index_2
end

A.data(:,2) = disp_v / 1000;

% make variable
time = A.data(:,1);
displacement =  A.data(:,2);
force =  A.data(:,3);
motor_speed = A.data(:,4);
motor_current = A.data(:,5);
strain =  A.data(:,6)/1000000;
stress = force/width/thickness;   

% Look at stress-strain at different strain pauses iE = 1:13
inds = [1492,11103; 
11226,20903;
21076,31003;
32266,41903;
42464,52603;
52809,63103;
63294,73504;
75454,88804;
88956,114854;
115036,126003;
127947,138203;
138397,148703;
148885,159273];

%% plot
clc;
close all;
clear tbl_cell;
for iE = 1:13
    ind_local = inds(iE,1):inds(iE,2);
    time_local = time(ind_local);
    time_local = time_local - time_local(1);
    stress_local = stress(ind_local);
    strain_local = strain(ind_local);
    
    figure; hold on;
    plot(time_local, strain_local);
    xlabel('Time (s)');
    ylabel('Strain');
    
    yyaxis right;
    plot(time_local, stress_local);
    ylabel('Stress (MPa)');
    set(gca,'fontsize',16);
    
    % estimate strain rate
    tbl = cell2table(cell(0,4));
    tbl.Properties.VariableNames = {'iE', 'itime', 'istress', 'strain_rate'};
    yyaxis left;
    
    nPts = 10; % nPts used to calculate strain rate, stress
    ind = 1;
    while ind+nPts-1 < length(time_local)
        if ind > 1000
            nPts = 1000;
        end
        ind_seg = ind: ind+nPts-1;
        stress_ii = mean(stress_local(ind_seg));
        strain_seg = strain_local(ind_seg);
        time_seg = time_local(ind_seg);
        time_ii = mean(time_seg);
        
        mdl = fitlm(time_seg, strain_seg);
        strain_rate = mdl.Coefficients.Estimate(2);    % rate = 1/s
        intercept = mdl.Coefficients.Estimate(1); 
        
        plot(time_seg, time_seg * strain_rate + intercept, '-', 'color', 'k', 'linewidth',3);
        tbl = [tbl; {iE, time_ii, stress_ii, strain_rate}];
        
        % increment ind
        ind = ind + 1.1 * nPts;
    end
    
    title(['iE=',num2str(iE)],'fontweight','normal');
    print(fullfile(save_dir, ['UM134 Mg_C3 iE=',num2str(iE),'.tiff']), '-dtiff');
    close;
    
    tbl_cell{iE} = tbl;
end


%% [2] tensile data
working_dir = 'E:\zhec umich Drive\2021-02-03 test relaxation using UM134 Mg_C3';

fileName = '2020-02-03 test relaxation.lvm';
delimiterIn = '\t';
headerlinesIn = 1;
A = importdata(fullfile(working_dir,fileName), delimiterIn, headerlinesIn);
display(A);

width = 3.45;       % UM134_Mg_C3, tested on 2021-01-29
thickness = 2.63;

% For this experiment, the displacement reading did not change well during test, so we need to manually fix
% For active loading part, the rate is 1 um/s. For other parts, the rate is 0 um/s  
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

% make variable
time = A.data(:,1);
displacement =  A.data(:,2);
force =  A.data(:,3);
motor_speed = A.data(:,4);
motor_current = A.data(:,5);
strain =  A.data(:,6)/1000000;
stress = force/width/thickness;   

% Look at stress-strain at different strain pauses iE = 1
inds(14,:) = [784, 13268]; 

%% plot
for iE = 14
    ind_local = inds(iE,1):inds(iE,2);
    time_local = time(ind_local);
    time_local = time_local - time_local(1);
    stress_local = stress(ind_local);
    strain_local = strain(ind_local);
    
    figure; hold on;
    plot(time_local, strain_local);
    xlabel('Time (s)');
    ylabel('Strain');
    
    yyaxis right;
    plot(time_local, stress_local);
    ylabel('Stress (MPa)');
    set(gca,'fontsize',16);
    
    % estimate strain rate
    tbl = cell2table(cell(0,4));
    tbl.Properties.VariableNames = {'iE', 'itime', 'istress', 'strain_rate'};
    yyaxis left;
    
    nPts = 10; % nPts used to calculate strain rate, stress
    ind = 1;
    while ind+nPts-1 < length(time_local)
        if ind > 1000
            nPts = 1000;
        end
        ind_seg = ind: ind+nPts-1;
        stress_ii = mean(stress_local(ind_seg));
        strain_seg = strain_local(ind_seg);
        time_seg = time_local(ind_seg);
        time_ii = mean(time_seg);
        
        mdl = fitlm(time_seg, strain_seg);
        strain_rate = mdl.Coefficients.Estimate(2);    % rate = 1/s
        intercept = mdl.Coefficients.Estimate(1); 
        
        plot(time_seg, time_seg * strain_rate + intercept, '-', 'color', 'k', 'linewidth',3);
        tbl = [tbl; {iE, time_ii, stress_ii, strain_rate}];
        
        % increment ind
        ind = ind + 1.1 * nPts;
    end
    
    title(['iE=',num2str(iE),' tensile'],'fontweight','normal');
    print(fullfile(save_dir, ['UM134 Mg_C3 tensile.tiff']), '-dtiff');
    close;
    
    tbl_cell{iE} = tbl;
end

%% Plot strain rate vs stress
close all;
for iE = 1:14
    legend_str{iE} = ['iE=',num2str(iE)];
    if iE==14
        legend_str{iE} = 'tensile';
    end
end

% plot strain rate of all points
figure; hold on;
colors = linspecer(14);
colororder(colors);
data_x = [];
data_y = [];
for iE = 1:14
    tbl = tbl_cell{iE};
    xx = tbl.istress;
    yy = tbl.strain_rate;
    plot(xx, yy, 'o', 'MarkerFaceColor',colors(iE,:), 'MarkerEdgeColor','k');
    data_x = [data_x; xx(:)];
    data_y = [data_y; yy(:)];
end
xlabel('stress (MPa)');
ylabel('strain rate, 1/s');
set(gca,'fontsize',16);
legend(legend_str,'location','eastoutside');

% plot strain rate of first few points
figure; hold on;
colors = linspecer(14);
colororder(colors);
data_x = [];
data_y = [];
for iE = 1:14
    tbl = tbl_cell{iE};
    xx = tbl.istress(1:1);
    yy = tbl.strain_rate(1:1);
    plot(xx, yy, 'o', 'MarkerFaceColor',colors(iE,:), 'MarkerEdgeColor','k');
    data_x = [data_x; xx(:)];
    data_y = [data_y; yy(:)];
end
xlabel('stress (MPa)');
ylabel('strain rate (1/s)');
set(gca,'fontsize',16);
legend(legend_str,'location','eastoutside');

% plot strain rate of last few points
figure; hold on;
colors = linspecer(14);
colororder(colors);
data_x = [];
data_y = [];
for iE = 1:14
    tbl = tbl_cell{iE};
    xx = tbl.istress(end-0:end);
    yy = tbl.strain_rate(end-0:end);
    plot(xx, yy, 'o', 'MarkerFaceColor',colors(iE,:), 'MarkerEdgeColor','k');
    data_x = [data_x; xx(:)];
    data_y = [data_y; yy(:)];
end
xlabel('stress (MPa)');
ylabel('strain rate (1/s)');
set(gca,'fontsize',16);
legend(legend_str,'location','eastoutside');


%% Log scale, using ABS of strain rate and stress
% close all;
% plot strain rate of first few points
figure; hold on;
colors = linspecer(14);
colororder(colors);
data_x = [];
data_y = [];
for iE = 1:14
    tbl = tbl_cell{iE};
    xx = log10(abs(tbl.istress(1:1)));
    yy = log10(abs(tbl.strain_rate(1:1)));
    plot(xx, yy, 'o', 'MarkerFaceColor',colors(iE,:), 'MarkerEdgeColor','k');
    data_x = [data_x; xx(:)];
    data_y = [data_y; yy(:)];
end
xlabel('log stress (MPa)');
ylabel('log strain rate (1/s)');
set(gca,'fontsize',16);
legend(legend_str,'location','eastoutside');
mdl = fitlm(data_x, data_y)
I = mdl.Coefficients.Estimate(1);
N = mdl.Coefficients.Estimate(2);
plot(min(data_x):max(data_x), [min(data_x):max(data_x)]*N + I, '-k');


% plot strain rate of last few points
figure; hold on;
colors = linspecer(14);
colororder(colors);
data_x = [];
data_y = [];
for iE = 1:14
    tbl = tbl_cell{iE};
    xx = log10(abs(tbl.istress(end-0:end)));
    yy = log10(abs(tbl.strain_rate(end-0:end)));
    plot(xx, yy, 'o', 'MarkerFaceColor',colors(iE,:), 'MarkerEdgeColor','k');
    data_x = [data_x; xx(:)];
    data_y = [data_y; yy(:)];
end
xlabel('log stress (MPa)');
ylabel('log strain rate (1/s)');
set(gca,'fontsize',16);
legend(legend_str,'location','eastoutside');
mdl = fitlm(data_x, data_y)
I = mdl.Coefficients.Estimate(1);
N = mdl.Coefficients.Estimate(2);
plot(min(data_x):0.1:max(data_x), [min(data_x):0.1:max(data_x)]*N + I, '-k');



