% Hank asked to compare experimental shear measured by AFM vs DIC

% Grain #1-9 in Hank's map have:
% ID = [80, 226, 185, 221, 83, 95, 189, 130, 176] in my EBSD data
clear;clc;
addChenFunction;
addpath('D:\p\m\CodeByDate');

% How to choose nDataPtsRange ========================================== 
% Assume DIC subsetSize = 2k+1, stepSize = s, filterSize = 2f+1
% A subset centered at ceil(k/s) will not be affected by the slipTraceLine 
% Due to averaging by filter, the non-affected subset center extended away from the slipTraceLine further by f data points
% Therefore, the non-affected subset center is ceil(k/s)+f data points away, or k+s*f pixels away. 
% In this study, nDataPtsRange > ceil(12/2)+2, so 10 seems fine.
nDataPtsRange = 10; % number of data points to cover on the horizontal line  

% Load orientation data
load('E:\Ti7Al_E1_insitu_tension\Analysis_r6c5\Ti7Al_E1_EbsdToSemForTraceAnalysis.mat', 'X','Y','boundaryTF', 'gID','gPhi1','gPhi','gPhi2','ID','eulerAligned'); 
load('D:\p\m\DIC_Analysis\setting_for_real_samples\Ti7Al_E1_setting.mat');

% Select area of interest (r6c5, r4c5, r7c5, ...) ============================
ir = 6;
ic = 5;

for iE = 1:7
    strainFile{iE} = matfile(['E:\Ti7Al_E1_insitu_tension\Analysis_r6c5\dic_s',num2str(iE),'.mat']);
end

%% ID: 226, grain-2, can select an iE 
iE = 7;

exx = strainFile{iE}.exx;
exy = strainFile{iE}.exy;
eyy = strainFile{iE}.eyy;
u = strainFile{iE}.u;
v = strainFile{iE}.v;
x = strainFile{iE}.x;
y = strainFile{iE}.y;

ID_current = 226;
ind = (gID==ID_current);
euler = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
[ssa, c_a, nss, ntwin, ssGroup] = define_SS('Ti','pyii'); % [plane normal; slip direction] 
[ss, c_a, ssa] = define_SS_cart('Ti','pyii');
m = angle2dcm(euler(1)/180*pi, euler(2)/180*pi, euler(3)/180*pi, 'zxz');

close all;
ss_to_compare = [4,5];  % selected slip systems to draw trace and compare  

myplot(x,y,exx,boundaryTF);title('\epsilon_x_x','fontweight','normal');
colormap(summer);
set(gca,'fontsize',18,'XTick',[0:1000:4000],'YTick',[0:1000:4000]);
xlabel('X (pixels)');
ylabel('Y (pixels)');

ir = 1060;
ic = 280;
rectangle('position',[2*ic,2*ir,2*340,2*340],'linewidth',3);
label_map_with_trace_for_Ti7Al_E1(x,y,ones(size(x))*ID_current,ID_current,ss_to_compare,'pyii',gca);    


%% can find the area of ~10um x 10um, which Hank analyzed by AFM
close all;
n_dp = 10 / (60/4096*2); % ~= 340 dpts, 680 pixels
ir = 1060;
ic = 280;
exx_local = exx(ir:ir+340,ic:ic+340);
% myplot(exx_local);
u = u(ir:ir+340,ic:ic+340);
v = v(ir:ir+340,ic:ic+340);
x = x(ir:ir+340,ic:ic+340); % x = x - min(x(:));
y = y(ir:ir+340,ic:ic+340); % y = y - min(y(:));
exx = exx(ir:ir+340,ic:ic+340);
exy = exy(ir:ir+340,ic:ic+340);
eyy = eyy(ir:ir+340,ic:ic+340);

myplot(x,y,exx);
set(gca,'fontsize',18);
title('\epsilon_x_x','fontweight','normal');
xlabel('X (pixels)');
ylabel('Y (pixels)');

myplot(x,y,u);
set(gca,'fontsize',18);
title('u','fontweight','normal');
xlabel('X (pixels)');
ylabel('Y (pixels)');

myplot(x,y,v)
set(gca,'fontsize',18);
title('v','fontweight','normal');
xlabel('X (pixels)');
ylabel('Y (pixels)');
%% Assume shear is due to activation of a single slip system
% (1): calculate 'shear'
% (2): back calculate step height ?
% Hank's tile size is 0.3um. 
% SEM-DIC is 60um, 4096 pixels, step size = 2, so 0.3 um covers
% 4096/2/60 * 0.3 = 34.13 * 0.3 = 10.24 data points, (happen to be similar to nDataPointsRange)  

% Method: 
% shear of ss = disp along Burgers / distance along planeNormal
% Assume using a virtual extensometer L = [dX, 0, 0]
% displacement S = [du,dv,dw] where dw = h(AFM) and is not known. 
% Assume S is due to activity of slip system: 
%   unit vector along slip direction b = [bx,by,bz]
%   slip plane normal n = [nx, ny, nz]
% Then, 
% (1) disp along Burgers S = du/bx = dv/by = dw/bz
% (2) distance along planeNormal = dot(L,n) = dX dot([1 0 0], n)

[nR,nC] = size(exx);
hws = 5;    % half window size

uL = zeros(nR,nC);
vL = zeros(nR,nC);
xL = zeros(nR,nC);
uR = zeros(nR,nC);
vR = zeros(nR,nC);
xR = zeros(nR,nC);
uL(:,1+hws:end-hws) = u(:,1:end-hws-hws);
vL(:,1+hws:end-hws) = v(:,1:end-hws-hws);
xL(:,1+hws:end-hws) = x(:,1:end-hws-hws);
uR(:,1+hws:end-hws) = u(:,1+hws+hws:end);
vR(:,1+hws:end-hws) = v(:,1+hws+hws:end);
xR(:,1+hws:end-hws) = x(:,1+hws+hws:end);

dU = uR - uL;   % unit: pixels
dV = vR - vL;
dX = xR - xL;   % For convenience, place virtual extensometer horizontally, so dY = 0, dZ = 0
dY = 0;

clear b n;
for ii = 4:5
    n{ii} = ss(1,:,ii) * m; % slip plane normal in sample coordinate
    b{ii} = ss(2,:,ii) * m; % burgers vector direction in sample coordinate
end
%% resolve du,dv onto ss#4 and #5
A = [1 1 0 0;
    0 0 1 1;
    b{4}(2), 0, -b{4}(1), 0;
    0, b{5}(2), 0, -b{5}(1)];

result = arrayfun(@(x,y) A\[x;y;0;0], dU, dV, 'UniformOutput', false);
du{4} = cellfun(@(x) x(1), result);
du{5} = cellfun(@(x) x(2), result);
dv{4} = cellfun(@(x) x(3), result);
dv{5} = cellfun(@(x) x(4), result);

myplot(dU);
set(gca,'fontsize',18);
xlabel('X (pixels)');
ylabel('Y (pixels)');
title('\Deltau','fontweight','normal');

myplot(x,y,du{4});
set(gca,'fontsize',18);
xlabel('X (pixels)');
ylabel('Y (pixels)');
title('\Deltau^4','fontweight','normal');

myplot(x,y,du{5});
set(gca,'fontsize',18);
xlabel('X (pixels)');
ylabel('Y (pixels)');
title('\Deltau^5','fontweight','normal');

myplot(dV);
set(gca,'fontsize',18);
xlabel('X (pixels)');
ylabel('Y (pixels)');
title('\Deltav','fontweight','normal');

myplot(x,y,dv{4});
set(gca,'fontsize',18);
xlabel('X (pixels)');
ylabel('Y (pixels)');
title('\Deltav^4','fontweight','normal');

myplot(x,y,dv{5});
set(gca,'fontsize',18);
xlabel('X (pixels)');
ylabel('Y (pixels)');
title('\Deltav^5','fontweight','normal');
%% repeat, but using du4 or du5 instead of du
for iss = 4:5
    % when we have 2 active systems, the system of equation is fully determined, so S1 and S2 should be the same  
    S1 = du{iss}/b{iss}(1);   % estimate by du, displacement along slip (Burgers vector) direction, unit: pixels
    S2 = dv{iss}/b{iss}(2);   % estimate by dv, displacement along slip (Burgers vector) direction

    NN = dX* dot([1 0 0],n{iss}); % distance along slip plane normal direction. A more precise expression is dot([dX,0,0], nn), as extensometer is along x-direction, so dY=0, dZ=0
    
    shear_1 = S1./NN;
    shear_2 = S2./NN;
    shear = (shear_1+shear_2)/2;  % shear, averaged of estimates by du and dv
    
    shearMap{iss} = shear_1;
    myplot(x,y,shearMap{iss});
    title(['\gamma^',num2str(iss)],'fontweight','normal');
    set(gca,'fontsize',18);
    xlabel('X (pixels)');
    ylabel('Y (pixels)');

end
shearT = shearMap{4} + shearMap{5};
myplot(x,y,shearT);
set(gca,'fontsize',18);
xlabel('X (pixels)');
ylabel('Y (pixels)');
title('\gamma','fontweight','normal');

shearT = my_pool(shearT,10,1);
xpool = my_pool(x,10,1);
ypool = my_pool(y,10,1);

myplot(xpool,ypool,shearT);
set(gca,'fontsize',18);
xlabel('X (pixels)');
ylabel('Y (pixels)');
title('\gamma','fontweight','normal');

%% If just consider one slip system
clear n b du dv;
for iss = 4:5 % select on system to analyze
    
    n{iss} = ss(1,:,iss) * m; % slip plane normal in sample coordinate
    b{iss} = ss(2,:,iss) * m; % burgers vector direction in sample coordinate
    
    S1 = dU/b{iss}(1);   % estimate by du, displacement along slip (Burgers vector) direction, unit: pixels
    S2 = dV/b{iss}(2);   % estimate by dv, displacement along slip (Burgers vector) direction
%     dw1 = S1*b{iss}(3) * (60/4096);  % sanity check: estimated by du, relative displacement in z-direction (dw or h, step height), unit: micron
%     dw2 = S2*b{iss}(3) * (60/4096);
    NN = dX* dot([1 0 0],n{iss}); % distance along slip plane normal direction. A more precise expression is dot([dX,0,0], nn), as extensometer is along x-direction, so dY=0, dZ=0
    
    shear_1 = S1./NN;
    shear_2 = S2./NN;
%     shear_1r = my_pool(shear_1,10,1);   % mean pool to make each tile = 10 data points, ~0.3 um
%     shear_2r = my_pool(shear_2,10,1);
    
    shear = (shear_1+shear_2)/2;  % shear, averaged of estimates by du and dv
%     shear_r = my_pool(shear,10,1);  % mean pool to make each tile = 10 data points, ~0.3 um
    
%     height = (dw1+dw2)/2;
%     height_r = my_pool(height,10,1);
    
    % Good: We can check if the estimates by du and dv are the same.
    % We can see that the difference is big, if the trace does not belong to this ss.
    shearMap{iss} = shear_1;
%     diffMap{iss} = abs(S1-S2);
    myplot(shear_1);
    myplot(shear_2);
    myplot(shear_1/2+shear_2/2);
    % myplot(diffMap{iss});
    
end












