
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

% (c) Number of realizations
Nr = 5;
figure;
for ii = 1:Nr
    [SIR,tVec] = ssir_model(tspan,X0);
    S_pred = SIR(:,1);
    I_pred = SIR(:,2);
    R_pred = SIR(:,3);
    C_pred = I_pred + R_pred;
    if ii == 1
        semilogy(tVec, C_pred, 'b',  tVec, R_pred, 'r', 1:length(C), C, 'bs', 1:length(R), R, 'r.');
        hold on;
    else
        semilogy(tVec, C_pred, 'b',  tVec, R_pred, 'r', 1:length(C), C, 'bs', 1:length(R), R, 'r.','HandleVisibility','off');
    end
    
end

title('sSIR model','fontweight','normal');
legend({'Cumulative infected (sSIR)','Deaths (sSIR)','Cumulative infected (data)','Deaths (data)'});
xlabel('Time (days since March 15');
ylabel('Number of individuals');
legend('location','best');
grid on;
set(gca,'fontsize',16);

