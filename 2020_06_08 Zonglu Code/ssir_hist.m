
close all;
clc;

covid_data = load('covid_data.mat').covid_data;   % load data
C = covid_data(:,1);    % Infected (cumulative)
R = covid_data(:,2);    % Removed (cumulative)

% initial condition
S0 = 500000;
I0 = 1;
R0 = 0;
X0 = [S0;I0;R0];
tspan = 100;

% (d) Number of realizations
Nr = 2000;
C = [];
for ii = 1:Nr
    [SIR,tVec] = ssir_model(tspan,X0);
    S_pred = SIR(:,1);
    I_pred = SIR(:,2);
    R_pred = SIR(:,3);
    C_pred = I_pred + R_pred;
    C(ii) = C_pred(end);
    disp(ii);
end

figure;
histogram(C, 37);

title('Histogram (total number of cases at day 100, 2000 realizations','fontweight','normal');
grid on;
set(gca,'fontsize',12);

% use cdf to estimate
figure;
[H, STATS] = cdfplot(C)
