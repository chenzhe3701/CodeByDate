covid_data = load('covid_data.mat').covid_data;   % load data
C = covid_data(:,1);    % Infected (cumulative) 
R = covid_data(:,2);    % Removed, i.e., Death, (cumulative)

days = 1:length(S);     % days

close all;
figure();
plot(days, C, 'b', days, R, 'r', 'linewidth', 3);     % plot and change linewidth to 3
legend({'Infected Individuals','Deaths'});     % add legend
xlabel('Time (days since March 15)');   % set axis label
ylabel('Number of individuals');
title('Covid-19 in Santa Barbara County','fontweight','normal');  % set title
legend('location','best');
grid on;    % turn grid on
set(gca,'fontsize',20); % set font size to 20

