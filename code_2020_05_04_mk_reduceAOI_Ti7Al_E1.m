% chenzhe, 2020-05-04
% based on the code for WE43_T6_C1 reduced AOI on 4/15/2020
% 
% The difference is: (1) the reduce AOI corresponds to an exact iRiC fov
% (2) the reduced AOI has higher resolution of DIC
% (3) There is no need to crop EBSD data, because of no gb adjustment;
% using EBSDtoSEM, we can directly read EBSD data, interp, and align.
%
% So, the main objective is to adjust the cpSEM field, and correct exy.

addChenFunction;
clear;
clc

%% load strain file, EBSD gb file, translation file during stitching
load('E:\Ti7Al_E1_insitu_tension\SEM Data\stitched_img\translations_searched_vertical_stop_0.mat','transX','transY');
load('E:\Ti7Al_E1_insitu_tension\Analysis_by_Matlab\Ti7Al_E1_EbsdToSemForTraceAnalysis.mat', 'ID','X','Y','boundaryTFB');
% % % load('E:\Ti7Al_E1_insitu_tension\Analysis_by_Matlab\Ti7Al_E1_EbsdToSemForTraceAnalysis.mat', 'ID','X','Y','boundaryTFB', ...
% % %     'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gNeighbors');
% % % target_EBSD_dir = uigetdir('E:\Ti7Al_E1_insitu_tension\Analysis_r7c5','dir to save EBSD for the reduced AOI');


% This is just plot the whole exx map for illustration
load('E:\Ti7Al_E1_insitu_tension\SEM Data\stitched_DIC\global_method_stitched_7.mat','exx'); 
myplot(X,Y,exx,boundaryTFB);

for iR = 1:8
   for iC = 1:8
       xt = transX(iR,iC);
       yt = transY(iR,iC);
       imrect(gca,[xt,yt,4096,4096]);
   end
end

%% So, select a FOV r6c5 is itself a good small AOI to analyze. Make new dic data
ir_AOI = 6;
ic_AOI = 5;

% source dir of high-res dic data for iRiC fov
dir_source = uigetdir(['E:\Ti7Al_E1_insitu_tension\Analysis_selected_area\r',num2str(ir_AOI),'c',num2str(ic_AOI)],'select source dir of the high-res DIC files');
% make new folder for the reduced AOI
target_dic_dir = uigetdir('E:\Ti7Al_E1_insitu_tension\Analysis_r6c5','select target dir to save updated DIC files');

% correct exy, correct sigma==-1 data
for iE = 0:7    
    file_source = ['Ti7Al_E1_S',num2str(iE),'_r',num2str(ir_AOI),'c',num2str(ic_AOI),'.mat'];
    load(fullfile(dir_source,file_source),'x','y','u','v','exx','exy','eyy','sigma');
    exy = -exy;
    exy_corrected = 1;
    exx(sigma==-1) = nan;
    exy(sigma==-1) = nan;
    eyy(sigma==-1) = nan;
    u(sigma==-1) = nan;
    v(sigma==-1) = nan;
    save([target_dic_dir,'\dic_s',num2str(iE),'.mat'],'x','y','u','v','exx','exy','eyy','sigma','exy_corrected','-v7.3');
end
%% copy and modify the setting file
[file_setting,dir_setting] = uigetfile('D:\p\m\DIC_Analysis\setting_for_real_samples\Ti7Al_E1_setting.mat','select the source setting file');
setting_file_source = [dir_setting,file_setting];
setting_file_new = [dir_setting,'\Ti7Al_E1_r',num2str(ir_AOI),'c',num2str(ic_AOI),'_setting.mat'];
copyfile(setting_file_source, setting_file_new);
load(setting_file_new,'cpSEM');
cpSEM(:,1) = cpSEM(:,1) - transX(1+ir_AOI, 1+ic_AOI);
cpSEM(:,2) = cpSEM(:,2) - transY(1+ir_AOI, 1+ic_AOI);
save(setting_file_new,'cpSEM','-append');

%% special case, if the EBSD data was processed to have gbAdjusted, to align with SEM-DIC, and want to read the adjusted EBSD data
% find inds in ID
% % % [~, indC_min] = min(abs(X(1,:)-(x(1,1) + transX(1+ir_AOI, 1+ic_AOI))));
% % % [~, indC_max] = min(abs(X(1,:)-(x(1,end) + transX(1+ir_AOI, 1+ic_AOI))));
% % % [~, indR_min] = min(abs(Y(:,1)-(y(1,1) + transY(1+ir_AOI, 1+ic_AOI))));
% % % [~, indR_max] = min(abs(Y(:,1)-(y(end,1) + transY(1+ir_AOI, 1+ic_AOI))));
% % % ID = ID(indR_min:indR_max, indC_min:indC_max);
% % % 
% % % %%
% % % save([target_EBSD_dir,'\EBSD_reduced_AOI.mat'], 'ID','gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gNeighbors');
    



