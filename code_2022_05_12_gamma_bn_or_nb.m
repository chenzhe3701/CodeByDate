%% gamma b n, Or, gamma n b
% Currently, my previous assumption is that F = I + du/dX = I + gamma * tensor_product( b(direction),  n(plane normal) ).  
% But try to see how different it is if n b.
% ==> comparing to experimental data, not so obvious. But on 2022-05-12, I
% am not so sure that my original assumption, bn is correct...

clear;
addChenFunction;
dicPath = uigetdir('D:\WE43_T6_C1\SEM Data\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('D:\p\m\DIC_Analysis\setting_for_real_samples\WE43_T6_C1_setting.mat','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1\Analysis_2021_09','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end
try
    load(fullfile(saveDataPath,[sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted']));
catch
    load(fullfile(saveDataPath,[sampleName,'_EbsdToSemForTraceAnalysis']));
end

gIDwithTrace = gID(~isnan(gExx));

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------

STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 6;

% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

neighbor_elim = 1;          % don't consider this ID as neighbor. For example, ID = 1 or 0 means bad region.
twinTF_text = 'twin';        % do you want to analyze twin? Use things like 'twin' or 'notwin'

data = matfile(fullfile(saveDataPath,[sampleName,'_3_twin_result_with_variant.mat']));

%% EDA
iE = 5;
vMap = data.variantMapCleanedCell;
vMap = vMap{iE};
strainFile = [dicPath,'\',f2,STOP{iE+B}]; disp(strainFile);
clear('exy_corrected');
load(strainFile,'exx','exy','eyy','sigma','exy_corrected');     % if 'exy_corrected' does not exist, this does not give error, rather, just warning. % ----------------------------------------------------------------------------------
if exist('exy_corrected','var')&&(1==exy_corrected)
    disp('================= exy already corrected ! ========================');
    exy_corrected = 1;
else
    disp('================= exy being corrected here ! =======================');
    exy = -exy;
    exy_corrected = 1;
end
% remove bad data points ----------------------------------------------------
exx(sigma==-1) = nan;
exy(sigma==-1) = nan;
eyy(sigma==-1) = nan;
qt_exx = quantile(exx(:),[0.0013,0.9987]); qt_exx(1)=min(-1,qt_exx(1)); qt_exx(2)=max(1,qt_exx(2));
qt_exy = quantile(exy(:),[0.0013,0.9987]); qt_exy(1)=min(-1,qt_exy(1)); qt_exy(2)=max(1,qt_exy(2));
qt_eyy = quantile(eyy(:),[0.0013,0.9987]); qt_eyy(1)=min(-1,qt_eyy(1)); qt_eyy(2)=max(1,qt_eyy(2));
ind_outlier = (exx<qt_exx(1))|(exx>qt_exx(2))|(exy<qt_exy(1))|(exy>qt_exy(2))|(eyy<qt_eyy(1))|(eyy>qt_eyy(2));
exx(ind_outlier) = nan;
exy(ind_outlier) = nan;
eyy(ind_outlier) = nan;

% myplot(exx, boundaryTFB);
% myplot(ID, boundaryTFB);

%% select grain
close all;
clc;

ID_current = 694;

% find index range of a small matrix containing the grain of interest,
% (and can choose to include some neighboring grains)
ind_local = ismember(ID, [ID_current]); %ismember(ID, [ID_current,ID_neighbor]);
indC_min = find(sum(ind_local, 1), 1, 'first');
indC_max = find(sum(ind_local, 1), 1, 'last');
indR_min = find(sum(ind_local, 2), 1, 'first');
indR_max = find(sum(ind_local, 2), 1, 'last');

exx_local = exx(indR_min:indR_max, indC_min:indC_max);  % strain of this region: grain + neighbor. Look at 'exx' strain, but can be changed later --------------------
exy_local = exy(indR_min:indR_max, indC_min:indC_max);
eyy_local = eyy(indR_min:indR_max, indC_min:indC_max);
boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
x_local = X(indR_min:indR_max, indC_min:indC_max);
y_local = Y(indR_min:indR_max, indC_min:indC_max);
ID_local = ID(indR_min:indR_max, indC_min:indC_max);
vMap_local = vMap(indR_min:indR_max, indC_min:indC_max);

myplot(vMap_local, boundaryTF_local); caxis([0 6]);
myplot(exx_local, boundaryTF_local);
myplot(exy_local, boundaryTF_local);
myplot(eyy_local, boundaryTF_local);

%%
[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
ss = crystal_to_cart_ss(ssa,c_a);

ind_euler = find(gID==ID_current);
euler = [gPhi1(ind_euler),gPhi(ind_euler),gPhi2(ind_euler)];
if (1==eulerAligned)
    g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
else
    g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
end
gamma = 0.1289; % twin shear for Mg
cPred = nan*zeros(nss,5);   % [iss, SF, exx, exy, eyy]
cPred_a = nan*zeros(nss,5);   % [iss, SF, exx, exy, eyy]
for iss = (nss+1):(nss+ntwin)   % for Mg
    %         disp('---');
    N(iss,:) = ss(1,:,iss) * g;
    M(iss,:) = ss(2,:,iss) * g;
    MN2{iss} = M(iss,:)'*N(iss,:);      %
    MN2{iss} = MN2{iss}(1:2,1:2);
    
    F3 = eye(3) + gamma * M(iss,:)' * N(iss,:);
    F = F3(1:2,1:2);
    % F = eye(2) + gamma*MN2{iss};
    epsilon = (F'*F-eye(2))/2;
    
    disp('----');
    disp((F3'*F3-eye(3))/2);
    % disp(epsilon);    

    cPred(iss,1) = iss;                                     % ss number
    cPred(iss,2) = N(iss,:) * stressTensor * M(iss,:)';     % Schmid factor
    cPred(iss,3:5) = [epsilon(1), epsilon(2), epsilon(4)];  % strain exx, exy, eyy.  Note that 'conjugated' twin system, i.e., 19 and 22, almost always show similar components!!!

    % Try alternative sequence
    FF3 = eye(3) + gamma * N(iss,:)' * M(iss,:);
    disp((FF3'*FF3-eye(3))/2);
    F = FF3(1:2,1:2);
    epsilon = (F'*F-eye(2))/2;
    cPred_a(iss,1) = iss;                                     % ss number
    cPred_a(iss,2) = N(iss,:) * stressTensor * M(iss,:)';     % Schmid factor
    cPred_a(iss,3:5) = [epsilon(1), epsilon(2), epsilon(4)];  % strain exx, exy, eyy.  Note that 'conjugated' twin system, i.e., 19 and 22, almost always show similar components!!!
end


% check strain of cluster #1, #4
c1 = 1;
c2 = 4;
cPred(19:end,:)
cPred_a(19:end,:)
inds = ismember(ID_local, ID_current) & ismember(vMap_local, c1);
cluster_1_mean = [nanmean(exx_local(inds)), nanmean(exy_local(inds)), nanmean(exy_local(inds))]
inds = ismember(ID_local, ID_current) & ismember(vMap_local, c2);
cluster_4_mean = [nanmean(exx_local(inds)), nanmean(exy_local(inds)), nanmean(exy_local(inds))]






