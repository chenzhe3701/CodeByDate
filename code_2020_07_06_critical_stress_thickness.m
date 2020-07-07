%% (1) setting for Ti-7Al, assume width = thickness
MAXLOAD = 4000; % max load can be applied
L = 19;   %   set sample length, mm
E = [117*1000, 2*1000, 1.273*1000];    % tangent modulus, MPa, lower
% E = [117*1000, 12.4*1000, 2*1000];    % tangent modulus, MPa, lower
Stress = [660, 770, 860];   % stress ranges, MPa
w0 = 2.2;     % minimum width of sample, mm

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
set(gca,'ylim',[0 1400]);
legend('buckling stress','max applicable stress')
xlabel('thickness,mm');ylabel('stress,MPa');title('Ti-7Al');

%% (2) setting for Mg WE43-T5
MAXLOAD = 3700; % max load can be applied
L = 28;   %   set sample length, mm
E = [45*1000, 5.55*1000, 0.707*1000, 0.374*1000];    % tangent modulus, MPa
Stress = [268, 283, 310, 320];   % stress ranges, MPa
w0 = 4.4;     % minimum width of sample, mm

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
xlabel('thickness,mm');ylabel('stress,MPa');title('WE43-T5');

%% (3) setting for Mg WE43-T6
MAXLOAD = 3700; % max load can be applied
L = 19;   %   set sample length, mm
E = [45*1000, 3.16*1000, 0.854*1000, 0.551*1000];    % tangent modulus, MPa
Stress = [160, 183, 211, 225];   % stress ranges, MPa
w0 = 4.5;     % minimum width of sample, mm

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

%% (4) setting for Mg WE43-T6
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