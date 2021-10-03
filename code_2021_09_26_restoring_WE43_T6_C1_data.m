% 2021-09-26
% Because of the disk crash and data loss
% (1) Reconstruct the 'struCell' for WE43_T6_C1 with correct 'cTrueTwin' field,
% (2) Then we can run 'remake_tNote_using_ref.m' to remake 'tNote'
% (3) Then we can run 'trace_analysis_3D_clusterToTwin_confirm', to quickly
% run analysis and get the repeatable result for confirming twin
% identificatoin.

%% 
% Load the published ground truth data of trueTwinMapCell
clear;
clc;
addChenFunction;

% ground truth data 'trueTwinMapCell'
load('D:\WE43_T6_C1\Analysis_2021_09\possibly useful data\WE43_T6_C1_trueTwinMapCell_published.mat','trueTwinMapCell');

% Load 'struCell' and 'twinMapCell' of initial identification
load('D:\WE43_T6_C1\Analysis_2021_09\20211001_0027_twinMaps.mat','struCell','twinMapCell');

% Load ID
saveDataPath = 'D:\WE43_T6_C1\Analysis_2021_09';
load(fullfile(saveDataPath,'WE43_T6_C1_EbsdToSemForTraceAnalysis_GbAdjusted.mat'), 'ID');

boundary = find_one_boundary_from_ID_matrix(ID);

% Comparing to the ground truth map, to construct 'tNote'
tNote = zeros(1,9);

%%
for iE = 2:5  % ==> iE
    % load clusterNumMapCleaned
    load(fullfile(saveDataPath,['WE43_T6_C1_s',num2str(iE),'_cluster_to_twin_result.mat']),'clusterNumMapCleaned');
    
    trueTwinMap = trueTwinMapCell{iE};  % 19:24
    twinMap = twinMapCell{iE};  % 1:6
    
    for iS = 1:length(struCell{iE})
        ID_current = struCell{iE}(iS).gID;   % ==> ID
        % iS = find(arrayfun(@(x) x.gID==179, struCell{iE}));
        disp(['ID = ',num2str(ID_current)]);
        % Initiate 'cTrueTwin'
        struCell{iE}(iS).cTrueTwin = struCell{iE}(iS).cActiveSS;
        
        nC = length(struCell{iE}(iS).cLabel);
        for iC = 1:nC    % ==> iC
            inds = ismember(ID,ID_current) & ismember(clusterNumMapCleaned, iC);
             
            tn_ground_truth = ismember(19:24, unique(trueTwinMap(inds)));   % ground trueth twin number
            % check the mode. If majority is 0, then the twin number might be due to clean up/ small inconsistency of maps?  
            ground_truth_mode = mode(trueTwinMap(inds));
            if ground_truth_mode == 0 
               tn_ground_truth = zeros(1,6); 
            end
            tn_current_label = ismember(1:6, unique(twinMap(inds)));        % currently labeled twin number
            
            if any(tn_ground_truth - tn_current_label)
                tNote = [tNote; ID_current, iE, iC, tn_ground_truth];
                fprintf('iE=%d, iD=%d, iC=%d \n', iE, ID_current, iC);
                disp(['size of tNote:' num2str(size(tNote))]);
            end
        end
    end
end




