% Analyze twinning activity for Mg4Al (sample S1 and C1)

clear;
addChenFunction;

%% Mg4Al_S1, iE = 0:6, iE=3 max compressive strain
output_dir = 'C:\Users\ZheChen\Desktop\Mg4Al Comparison';

% looks like have to include this part to read the sample name.
path_setting = 'D:\p\m\DIC_Analysis\setting_for_real_samples\';
file_setting = 'Mg4Al_S1_setting.mat';
load_settings([path_setting,file_setting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor','strainPauses');

data_dir = 'E:\Mg4Al_S1_insitu\Analysis_by_Matlab\';
% load([data_dir,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');
data_file = matfile([data_dir,sampleName,'_EbsdToSemForTraceAnalysis']);
load([data_dir,sampleName,'_EbsdToSemForTraceAnalysis'],'gID','gPhi1','gPhi','gPhi2');

twin_variant_data = 'D:\p\m\DIC_Analysis\temp_results\20200325_0506_new_variant_map_Mg4Al_S1.mat';
load(twin_variant_data, 'struCell','variantMapCleanedCell')

um_per_dp = 5*60/4096;

% Table_1, for grain summaries.
variableNames = {'ID','twin_SF','activeTwin_SF','gDia','twinnedTF','tAF','euler'};
T = cell2table(cell(0,length(variableNames)));
T.Properties.VariableNames = variableNames;

% Table_2, summarize things for each of the 6 possible variant, for all grains
variableNames2 = {'ID','iTwin','variant_SF','vActiveTF','vPct'};
T2 = cell2table(cell(0,length(variableNames2)));
T2.Properties.VariableNames = variableNames2;

for iE = 3
    for iS = 1:length(struCell{iE})
        ID_current = struCell{iE}(iS).gID;
        gd = sqrt(4*struCell{iE}(iS).gVol/pi) * um_per_dp;
        vol_t = sum(struCell{iE}(iS).tVol);
        vol_g = struCell{iE}(iS).gVol;
        
        ind = find(ID_current==gID);
        euler_current = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
        [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler_current, [0 0 0], [0 0 0], stressTensor, 'Mg', 'twin');
        twin_SF = max(abs_schmid_factor(19:24,2));
        
        activeTS = logical(sum(struCell{iE}(iS).cTrueTwin,1));
        tSFs = struCell{iE}(iS).tSF;
        
        if any(activeTS)
            activeTwin_SF = max(struCell{iE}(iS).tSF(activeTS));
            T = [T; {ID_current, twin_SF, activeTwin_SF, gd, true, vol_t/vol_g, euler_current}];
        else
            activeTwin_SF = nan;
            T = [T; {ID_current, twin_SF, activeTwin_SF, gd, false, vol_t/vol_g, euler_current}];
        end
        
        for iTwin=1:length(activeTS)
            if activeTS(iTwin)>0
                vActiveTF = true;
                T2 = [T2; {ID_current, iTwin, tSFs(iTwin), vActiveTF, struCell{iE}(iS).tVol(iTwin)/vol_g}];
            else
                vActiveTF = false;
                T2 = [T2; {ID_current, iTwin, tSFs(iTwin), vActiveTF, struCell{iE}(iS).tVol(iTwin)/vol_g}];
            end
        end
        
    end
end

ind = (T.twinnedTF==1);
eulers_twinned = T.euler(ind,:);
eulers_not_twinned = T.euler(~ind,:);

plot_on_IPF(eulers_twinned, [0 0 0], [0 0 0], [1 0 0], [19:24], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
title('Twinned');
plot_on_IPF(eulers_not_twinned, [0 0 0], [0 0 0], [1 0 0], [19:24], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
title('Not twinned');

%% Mg4Al_C1, iE = 0:7, iE=4 max strain

data_dir = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD\analysis\';

% Table_1, for grain summaries.
variableNames = {'ID','twin_SF','twinnedTF','tAF','euler'};
T = cell2table(cell(0,length(variableNames)));
T.Properties.VariableNames = variableNames;

% Table_2, summarize things for each of the 6 possible variant, for all grains
variableNames2 = {'ID','iTwin','variant_SF','vActiveTF'};
T2 = cell2table(cell(0,length(variableNames2)));
T2.Properties.VariableNames = variableNames2;

load([data_dir,'variant_maps'],'variantMap');
% load([data_dir,'geotrans_and_id_link'],'tbl'); % table width = 8, iE=0 at col 1 --> iE at col iE+1

iE = 4;
% with these data overwrite, no need to check the linkage. But valid grain
% IDs are between 1 and 1000
load([data_dir,'data_with_ID_overwrite_iE_',num2str(iE),'.mat'], 'ID','gID','gPhi1','gPhi','gPhi2');
[gPhi1, gPhi, gPhi2] = align_euler_to_sample(gPhi1, gPhi, gPhi2, 'custom', 90, 180, 0); % setting-1 align

inds = (gID>0) & (gID<1000);
valid_gID = gID(inds);

variant_map = variantMap{iE};
for ii = 1:length(valid_gID)
    ID_current = valid_gID(ii);
    
    vol_t = sum((ID(:)==ID_current)&(variant_map(:)>0));
    vol_g = sum(ID(:)==ID_current);        
        
    ind = find(ID_current==gID);
    euler_current = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
    
    [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler_current, [0 0 0], [0 0 0], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
    tSFs = abs_schmid_factor(19:24,2);
    twin_SF = max(tSFs);
    
    active_variants = unique(variant_map(ID==ID_current));
    active_variants(active_variants == 0) = [];
    activeTS = false(1,6);  % logical representing if variant is active
    activeTS(active_variants) = 1;
    
    if any(activeTS)
        T = [T; {ID_current, twin_SF, true, vol_t/vol_g, euler_current}];
    else
        T = [T; {ID_current, twin_SF, false, vol_t/vol_g, euler_current}];
    end
    
    for iTwin=1:length(activeTS)
        if activeTS(iTwin)>0
            vActiveTF = true;
            T2 = [T2; {ID_current, iTwin, tSFs(iTwin), vActiveTF}];
        else
            vActiveTF = false;
            T2 = [T2; {ID_current, iTwin, tSFs(iTwin), vActiveTF}];
        end
    end
        
end

ind = (T.twinnedTF==1);
eulers_twinned = T.euler(ind,:);
eulers_not_twinned = T.euler(~ind,:);

plot_on_IPF(eulers_twinned, [0 0 0], [0 0 0], [1 0 0], [19:24], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
title('Twinned');
plot_on_IPF(eulers_not_twinned, [0 0 0], [0 0 0], [1 0 0], [19:24], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
title('Not twinned');




% variant area fraction
ind = (T2.iE==iE)&(T2.vActiveTF==1);
t2 = T2(ind,:);
edges = 0:0.05:0.5;
vx = t2.variant_SF;                   % ---------------------------> Select which SF here
vy = t2.vPct;                      % ---------------------------> Select strain here
ylabel_str = 'Twin Variant Area Fraction in Grain';    % ---------------------------> change name

figure;disableDefaultInteractivity(gca);
plot(vx, vy, '.');
xlabel('Schmid Factor'); ylabel(ylabel_str);
% set(gca,'ylim',[0 0.1])
title(['iE=',num2str(iE),', twinned grains'],'fontweight','normal');

gv = discretize(vx, edges);
nGroups = length(edges)-1;
clear labels;
for ii = 1:length(edges)-1
    labels{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
end
figure;disableDefaultInteractivity(gca);
boxplot([vy; nan*ones(nGroups, 1)], [gv; (1:nGroups)'],'notch','on');
xlabel('Twin Variant Schmid Factor'); ylabel(ylabel_str)
set(gca,'xticklabels',labels,'xticklabelrotation',45,'fontsize',15);
title(['iE=',num2str(iE),', twinned grains'],'fontweight','normal');
