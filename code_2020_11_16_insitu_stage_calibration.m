% The microstrain reading from NI-DAQ, converted from analog output of the
% vishay model p3 strain gage reader, needs to be calibrated.

working_dir = 'E:\zhec umich Drive\2020-11-13 Mg4Al_C2 insiu curve';
cd(working_dir);

%% read data
fileName = '2020-11-14 insitu stage calibration.xlsx';

T = readtable(fullfile(working_dir, fileName));

data = T.Variables;
micro_strain = data(:,1);
voltage = data(:,2)/10000;

%%
close all;
figure;
plot(voltage, micro_strain, '-o');
xlabel('voltage (V)');
ylabel('micro strain');
set(gca,'fontsize',18);

% micro_strain = k * voltage + b = k * (voltage + b/k);
mdl = fitlm(voltage, micro_strain);
k = mdl.Coefficients.Estimate(2);
b = mdl.Coefficients.Estimate(1);
b_over_k = b/k;

str = sprintf('micro strain = %.0f x (voltage %.4f)', k, b_over_k); 
legend(str,'location','northoutside');

% result: micro_strain = 24382 x (voltage - 1.2606)