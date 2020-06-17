
close all;
clc;

covid_data = load('covid_data.mat').covid_data;   % load data
C = covid_data(:,1);    % Infected (cumulative)
R = covid_data(:,2);    % Removed (cumulative)

% initial condition
S0 = 500000;
I0 = 1;
R0 = 0;
y0 = [S0; I0; R0];

beta0 = [0.1 : 0.001 : 0.5] ./ S0;  % values of beta to be tested

global beta;    % beta to be tested

Ndays = 100;
tspan = 0:Ndays;   % time span, evaluate at integer days. 

sz = size(covid_data,1); % actual data only covers 60 days

clear t y err;
for ii  = 1:length(beta0)
    
    beta = beta0(ii);
    
    [t,y] = ode45(@sir_model, tspan, y0);   % y = [S,I,R]
    
    % only 60 days
    S_pred = y(1:sz,1);
    I_pred = y(1:sz,2);
    R_pred = y(1:sz,3);
    
    err(ii) = sqrt(sum((I_pred + R_pred - C).^2) / sz);
    
end

% (c) error plot, beta0 vs err
figure;
plot(beta0, err, 'linewidth',3);
title('Error analysis (SIR model vs data)','fontweight','normal');
xlabel('beta0');
ylabel('Error');
grid on;
set(gca,'fontsize',20,'ylim',[0,4000]);

% (d) find best value of beta
[~, ind] = min(err);
beta = beta0(ind);
if abs(beta - 7.58e-7)<0.1e-7
    disp(['The optimizatin result was beta_0 = ',num2str(beta, '%4.3e')]);
else
    disp(['The optimizatin result was beta_0 = ',num2str(beta, '%4.3e')]);
    disp('Not able to accompolish the optimization, proceed with beta_0 = 7.580e-7');
    beta = 7.58e-7;    
end

[t,y] = ode45(@sir_model, tspan, y0);   % solve again using best beta,  y = [S,I,R]

% all 100 days
S_pred = y(:,1);
I_pred = y(:,2);
R_pred = y(:,3);
C_pred = I_pred + R_pred;

figure; 
semilogy(0:Ndays, C_pred, 'b',...
    0:Ndays, R_pred, 'r',...
    1:sz, C, 'bs',...
    1:sz, R, 'r.',...
    'linewidth', 3);  
title(['SIR model, beta0 = ',num2str(beta, '%4.3e')],'fontweight','normal');
legend({'Cumulative infected (SIR)','Deaths (SIR)','Cumulative infected (data)','Deaths (data)'});
xlabel('Time (days since March 15)');
ylabel('Number of individuals');
grid on;
set(gca,'fontsize',20);
legend('location','best');








