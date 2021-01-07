
%% test loading using labview and new insitu stage

%% import data. Can change to #1, #2, #3
filename = 'data\tensile aluminum load to 500N #1.lvm';
delimiterIn = '\t';
headerlinesIn = 1;
A = importdata(filename, delimiterIn, headerlinesIn);
display(A.colheaders)

[~, ind] = max(A.data(:,3));
A.data = A.data(1:ind,:);

%% make variable
time = A.data(:,1);
displacement = A.data(:,2);
force = A.data(:,3);
motor_speed = A.data(:,4);
motor_current = A.data(:,5);
strain = A.data(:,6)/1000000;
stress = force/1.75/2.95;   % pure Mg sample, cross section
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
plot(time, motor_speed, '-r', 'linewidth', 3);
xlabel('time (sec)');
ylabel('motor speed indication (V)');
set(gca,'fontsize',16);

figure;
plot(time, motor_current, '-r', 'linewidth', 3);
xlabel('time (sec)');
ylabel('motor current indication (V)');
set(gca,'fontsize',16);

figure;
plot(time, strain, '-r', 'linewidth', 3);
xlabel('time (sec)');
ylabel('strain');
set(gca,'fontsize',16);

%% stress-strain relationship
close all;
figure;
plot(strain, stress, '-r', 'linewidth', 3);
set(gca,'ylim',[-20 100]);
xlabel('strain');
ylabel('stress (MPa)');
set(gca,'fontsize',16);

nPts = round(ind * 0.8);
md = fitlm(strain(1:nPts), stress(1:nPts));
b = md.Coefficients.Estimate(1);
k = md.Coefficients.Estimate(2);
drawline('Position',[0,b; strain(nPts), strain(nPts)*k+b]);
text(strain(nPts)/2, (strain(nPts)*k+b -20)/2, ['E = ',num2str(k/1000, '%.1f'), ' GPa'], 'fontsize', 16);

%% load vs motor_current
figure;
plot(force, motor_current, '-r', 'linewidth', 3);
xlabel('force (N)');
ylabel('motor current indication (V)');
set(gca,'fontsize',16);



