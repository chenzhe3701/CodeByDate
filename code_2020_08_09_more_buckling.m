clear;

%% Assume Euler buckling:
% The critical load is: P_cr = 4*pi^2*E*I/L^2
% Therefore, sigma_cr = P_cr/A = E * I/A * 4*pi^2/L^2 = E * (w*t^3/12)/(wt) *(4*pi^2/L^2)
% = (E*pi^2)/3 * (t^2/L^2)

MAXLOAD = 1000; % max load can be applied
L = 20;   %   set sample length, mm
w0 = 2;     % minimum width of sample, mm

%% (1) PureMg data from Aeirel
% import stress strain data
% format: [time(sec), disp(mm), force(N), stress(MPa), strain(%)]
data = readmatrix('data\Pure_Mg_Tensile_Data.xlsx');
cross_section = 28.218; % mm^2

stress_initial = data(:,4);    % 
stress = data(:,3)/cross_section;  %

strain = data(:,5);    % Percent
strain = strain / 100;

[~, ind_end] = min(abs(strain-0.08)) % analyze up to 0.08 strain

% illustrate stress, strain, find useful portion
close all;
figure; hold on;
% plot(strain, stress_initial,'r.-','linewidth',4);
plot(strain(1:ind_end), stress(1:ind_end),'b.-');
xlabel('Strain');
ylabel('Stress (MPa)');
title('Pure Mg');
set(gca,'fontsize',18);



%% Euler buckling:
% The critical load is: P_cr = 4*pi^2*E*I/L^2
% Therefore, sigma_cr = P_cr/A = E * I/A * 4*pi^2/L^2 = E * (w*t^3/12)/(wt) *(4*pi^2/L^2)
% = (E*pi^2)/3 * (t^2/L^2)

clear E Stress t;

npts = 500;
jj = 1;
for ii = 1:round(ind_end/npts)
    stress_part = stress(1+(ii-1)*npts : npts*ii);
    strain_part = strain(1+(ii-1)*npts : npts*ii);
    md = fitlm(strain_part, stress_part);
    if md.Coefficients.Estimate(2)>0
        E(jj) = md.Coefficients.Estimate(2);    % tangent modulus, MPa, lower
        Stress(jj) = mean(stress_part);         % stress ranges, MPa
        jj = jj + 1;
    end
end

for ii = 1:length(E)
    t(ii) = sqrt(Stress(ii)*3*L^2/E(ii)/pi^2);
end

% figure;
% plot(t,Stress)

w = w0*ones(size(t));
w(w<t) = t(w<t);   % we want to make w>=thick
sigmaMax = MAXLOAD./t./w;
figure;
plot(t,Stress,'.b-',t,sigmaMax,'.r');
set(gca,'ylim',[0 300], 'xlim', [0,10]);
legend('buckling stress','max applicable stress')
xlabel('thickness,mm');ylabel('stress,MPa');title('Pure Mg');
title('Pure Mg');
set(gca,'fontsize',18);

%% (2) Mg4Al_S1 compresion data
% import stress strain data
data = load('data\Mg4Al_S1_compression_curve_data.mat');
stress = -data.stress_r;  %
% strain = -data.exx_interped_by_dic_r;
strain = -data.strain_r;
[~, ind_end] = min(abs(strain-0.025)); % analyze up to 0.08 strain. No need for WE43-T6, but doesn't matter. 

% illustrate stress, strain, find useful portion
figure; hold on;
plot(strain(1:ind_end), stress(1:ind_end),'b.-');
xlabel('Strain');
ylabel('Stress (MPa)');
title('Mg4Al S1');
set(gca,'fontsize',18);


%% Euler buckling:
% The critical load is: P_cr = 4*pi^2*E*I/L^2
% Therefore, sigma_cr = P_cr/A = E * I/A * 4*pi^2/L^2 = E * (w*t^3/12)/(wt) *(4*pi^2/L^2)
% = (E*pi^2)/3 * (t^2/L^2)

clear E Stress t;

npts = 40;
jj = 1;
for ii = 1:floor(ind_end/npts)
    stress_part = stress(1+(ii-1)*npts : npts*ii);
    strain_part = strain(1+(ii-1)*npts : npts*ii);
    md = fitlm(strain_part, stress_part);
    if md.Coefficients.Estimate(2)>0
        E(jj) = md.Coefficients.Estimate(2);    % tangent modulus, MPa, lower
        Stress(jj) = mean(stress_part);         % stress ranges, MPa
        jj = jj + 1;
    end
end

for ii = 1:length(E)
    t(ii) = sqrt(Stress(ii)*3*L^2/E(ii)/pi^2);
end

% figure;
% plot(t,Stress)

w = w0*ones(size(t));
w(w<t) = t(w<t);   % we want to make w>=thick
sigmaMax = MAXLOAD./t./w;
figure;
plot(t,Stress,'.b-',t,sigmaMax,'.r--');
set(gca,'ylim',[0 150], 'xlim', [0,7]);
legend('buckling stress','max applicable stress')
xlabel('thickness,mm');ylabel('stress,MPa');title('WE43-T6');
title('Mg4Al');
set(gca,'fontsize',18);


%% (3) WE43-T6 compresion data
% import stress strain data
data = load('data\WE43_T6_C1_compression_curve_data.mat');
stress = -data.stress_r;  %
% strain = -data.exx_interped_by_dic_r;
strain = -data.strain_r;

% illustrate stress, strain, find useful portion
figure; hold on;
plot(strain, stress,'b.-');
xlabel('Strain');
ylabel('Stress (MPa)');
title('WE43-T6');
set(gca,'fontsize',18);
ind_end = length(strain); % analyze up to 0.08 strain. No need for WE43-T6, but doesn't matter. 

%% Euler buckling:
% The critical load is: P_cr = 4*pi^2*E*I/L^2
% Therefore, sigma_cr = P_cr/A = E * I/A * 4*pi^2/L^2 = E * (w*t^3/12)/(wt) *(4*pi^2/L^2)
% = (E*pi^2)/3 * (t^2/L^2)

clear E Stress t;

npts = 50;
jj = 1;
for ii = 1:floor(ind_end/npts)
    stress_part = stress(1+(ii-1)*npts : npts*ii);
    strain_part = strain(1+(ii-1)*npts : npts*ii);
    md = fitlm(strain_part, stress_part);
    if md.Coefficients.Estimate(2)>0
        E(jj) = md.Coefficients.Estimate(2);    % tangent modulus, MPa, lower
        Stress(jj) = mean(stress_part);         % stress ranges, MPa
        jj = jj + 1;
    end
end

for ii = 1:length(E)
    t(ii) = sqrt(Stress(ii)*3*L^2/E(ii)/pi^2);
end

% figure;
% plot(t,Stress)

w = w0*ones(size(t));
w(w<t) = t(w<t);   % we want to make w>=thick
sigmaMax = MAXLOAD./t./w;
figure;
plot(t,Stress,'.b-',t,sigmaMax,'.r--');
set(gca,'ylim',[0 600], 'xlim', [0,7]);
legend('buckling stress','max applicable stress')
xlabel('thickness,mm');ylabel('stress,MPa');title('WE43-T6');
title('WE43-T6');
set(gca,'fontsize',18);

%%
if 0
    %% (3.1) Set interval the same as original calculation for comparison to double check
    
    clear E Stress t;
    
    npts = 200;
    jj = 1;
    for ii = 1:floor(ind_end/npts)
        stress_part = stress(1+(ii-1)*npts : npts*ii);
        strain_part = strain(1+(ii-1)*npts : npts*ii);
        md = fitlm(strain_part, stress_part);
        if md.Coefficients.Estimate(2)>0
            E(jj) = md.Coefficients.Estimate(2);    % tangent modulus, MPa, lower
            Stress(jj) = mean(stress_part);         % stress ranges, MPa
            jj = jj + 1;
        end
    end
    
    MAXLOAD = 3700; % max load can be applied
    L = 19;   %   set sample length, mm
    w0 = 4.5;     % minimum width of sample, mm
    
    for ii = 1:length(E)
        t(ii) = sqrt(Stress(ii)*3*L^2/E(ii)/pi^2);
    end
    
    % figure;
    % plot(t,Stress)
    
    w = w0*ones(size(t));
    w(w<t) = t(w<t);   % we want to make w>=thick
    sigmaMax = MAXLOAD./t./w;
    figure;
    plot(t,Stress,'.b--',t,sigmaMax,'.r--');
    set(gca,'ylim',[0 600], 'xlim', [0,7]);
    legend('buckling stress','max applicable stress')
    xlabel('thickness,mm');ylabel('stress,MPa');title('WE43-T6');
    
    %% (3.2) The original one used for WE43-T6 for comparison
    MAXLOAD = 3700; % max load can be applied
    L = 18;   %   set sample length, mm
    E = [45*1000, 20*1000, 7*1000, 0.551*1000];    % tangent modulus, MPa
    Stress = [160, 180, 222, 225];   % stress ranges, MPa
    w0 = 3.5;     % minimum width of sample, mm
    
    % plot
    clear sigma;
    clear t;
    t_min_required = 0;     % min_required thickness, initialize = 0
    for ii=1:length(E)
        % target stress
        if ii==1
            sigma{ii} = linspace(0,Stress(ii),50);
        else
            sigma{ii} = linspace(Stress(ii-1),Stress(ii),50);
        end
        t{ii} = sqrt(sigma{ii}*3*L^2/E(ii)/pi^2);   % buckling thickness
        
        ind = (t{ii} >= t_min_required);
        t{ii} = t{ii}(ind);
        sigma{ii} = sigma{ii}(ind);
        t_min_required = max(cell2mat(t));  % min required thickness, after this step
    end
    
    t = cell2mat(t);
    sigma = cell2mat(sigma);
    
    t2 = 0.01:0.01:max(t);
    w = w0*ones(size(t2));
    w(w<t2) = t2(w<t2);   % we want to make w>=thick
    sigmaMax = MAXLOAD./t2./w;
    figure;
    plot(t,sigma,'-o',t2,sigmaMax);
    set(gca,'ylim',[0 600]);
    legend('buckling stress','max applicable stress')
    xlabel('thickness,mm');ylabel('stress,MPa');title('WE43-T6');
    
end