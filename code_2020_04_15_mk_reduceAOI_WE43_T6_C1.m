% chenzhe, 2020-04-15
% create a small AOI from WE4_T6_C1, to check codes/algorithms quickly
clear;
clc

%% load strain file, EBSD gb file, translation file during stitching
load('D:\WE43_T6_C1\SEM Data\stitched_DIC\_5.mat','exx');
load('D:\WE43_T6_C1\SEM Data\stitched_img\translations_searched_vertical_stop_0.mat','transX','transY');
load('D:\WE43_T6_C1\Analysis_by_Matlab_after_realign\WE43_T6_C1_EbsdToSemForTraceAnalysis_GbAdjusted_v73.mat', 'ID','X','Y','boundaryTFB', ...
    'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gNeighbors');
target_EBSD_dir = uigetdir('D:\WE43_T6_C1\Analysis_by_Matlab_reducedAOI','dir for the reduced AOI EBSD data');
% plot
myplot(X,Y,exx,boundaryTFB);

for iR = 1:11
   for iC = 1:19
       xt = transX(iR,iC);
       yt = transY(iR,iC);
       imrect(gca,[xt,yt,4096,4096]);
   end
end

%% make new folder for the reduced AOI
dir_target = uigetdir('D:\WE43_T6_C1\SEM Data','select dir where you want to create a folder for data of the selected AOI');
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(dir_target,'DIC reduced AOI')
if SUCCESS
    dir_target = [dir_target,'\DIC reduced AOI'];
else
    error('dir not succesfully made');
end

%% So, looks like the FOV r4c18 is itself a good small AOI to analyze. Make new dic data
ir_AOI = 4;
ic_AOI = 18;
for iE = 0:7
    dir_source = ['D:\WE43_T6_C1\SEM Data\byFov\r',num2str(ir_AOI),'c',num2str(ic_AOI)];
    file_source = ['WE43_T6_C1_s',num2str(iE),'_r',num2str(ir_AOI),'c',num2str(ic_AOI),'.mat'];
    load(fullfile(dir_source,file_source),'x','y','u','v','exx','exy','eyy','sigma');
    exy = -exy;
    exy_corrected = 1;
    save([dir_target,'\_',num2str(iE),'.mat'],'x','y','u','v','exx','exy','eyy','sigma','exy_corrected');
end
%% copy and modify the setting file
dir_setting = uigetdir('D:\p\m\DIC_Analysis\setting_for_real_samples\','select dir for the setting files');
setting_file_source = [dir_setting,'\WE43_T6_C1_setting.mat'];
setting_file_new = [dir_setting,'\WE43_T6_C1_reducedAOI_setting.mat'];
copyfile(setting_file_source, setting_file_new);
load(setting_file_new,'cpSEM');
cpSEM(:,1) = cpSEM(:,1) - transX(1+ir_AOI, 1+ic_AOI);
cpSEM(:,2) = cpSEM(:,2) - transY(1+ir_AOI, 1+ic_AOI);
save(setting_file_new,'cpSEM','-append');

%% special case, if the EBSD data was processed to have gbAdjusted, to align with SEM-DIC, and want to read the adjusted EBSD data
% find inds in ID
[~, indC_min] = min(abs(X(1,:)-(x(1,1) + transX(1+ir_AOI, 1+ic_AOI))));
[~, indC_max] = min(abs(X(1,:)-(x(1,end) + transX(1+ir_AOI, 1+ic_AOI))));
[~, indR_min] = min(abs(Y(:,1)-(y(1,1) + transY(1+ir_AOI, 1+ic_AOI))));
[~, indR_max] = min(abs(Y(:,1)-(y(end,1) + transY(1+ir_AOI, 1+ic_AOI))));
ID = ID(indR_min:indR_max, indC_min:indC_max);
save([target_EBSD_dir,'\EBSD_reduced_AOI.mat'], 'ID','gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gNeighbors');
    



