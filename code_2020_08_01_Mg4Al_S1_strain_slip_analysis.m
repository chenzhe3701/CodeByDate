%% try to modify and organize this code for paper
% (1) Identify slip traces (mainly basal + some prismatic)

clear;
clc;
addChenFunction;

% setting file
fileSetting = 'Mg4Al_S1_setting.mat';
pathSetting = 'D:\p\m\DIC_Analysis\setting_for_real_samples';
sampleName = [];    % such as 'Ti7Al_#B6'
sampleMaterial = [];  % such as 'Ti' or 'Mg'
stressTensor = [];
load_settings(fullfile(pathSetting,fileSetting),'sampleName','sampleMaterial','stressTensor','strainPauses');

% path for saved data
saveDataPath = 'E:\Mg4Al_S1_insitu\Analysis_by_Matlab';
saveDataFile = fullfile(saveDataPath,'Mg4Al_S1_EbsdToSemForTraceAnalysis.mat');
cd(saveDataPath);

dataFile = matfile(saveDataFile)
load(saveDataFile,'ID','X','Y','boundaryTF','boundaryTFB','gDiameter','gID','gNeighbors','gNNeighbors','gPhi1','gPhi','gPhi2','eulerAligned');

% strain data
dicDataPath = 'E:\Mg4Al_S1_insitu\SEM Data\stitched_DIC';
for iE = 1:6
    dicFile{iE} = fullfile(dicDataPath,['_',num2str(iE),'_v73.mat']);
    dicData{iE} = matfile(dicFile{iE});
end

% struCell
variant_file = 'D:\p\m\DIC_Analysis\temp_results\20200325_0506_new_variant_map_Mg4Al_S1.mat';
load(variant_file, 'struCell', 'trueTwinMapCell');

output_dir  = 'E:\zhec umich Drive\0_temp_output\Mg4Al_S1 analysis\for paper';
mkdir(output_dir);
%% calculate effective strain
for iE = 1:6
    eEff{iE} = calculate_effective_strain(dicData{iE}.exx, dicData{iE}.exy, dicData{iE}.eyy);
end

%% Exploritory analysis. Find how many grains in the AOI. Plot grain size distribution in AOI.
disp(['min(gID) = ', num2str(min(gID))]);    % smallest gID
disp(['max(gID) = ', num2str(max(gID))]);    % max gID
disp(['numel(gID) = ', num2str(numel(gID))]);    % number of grains
gID_list = unique(ID(:));   % grains that exist in the AOI
disp(['numel(gID_list) = ', num2str(numel(gID_list))]); % number of grains exist in the AOI

%% struCell has only 177 grains, gID_list has 180. find the difference
ID_in_struCell = [];
for ii = 1:length(struCell{2})
    ID_in_struCell = [ID_in_struCell, struCell{2}(ii).gID];
end
ID_diff = gID_list(~ismember(gID_list, ID_in_struCell));
myplot(ismember(ID,ID_diff));
gDiameter(ismember(gID,ID_diff))
% -> Note: the result shows that the 3 grains are edge grains on EBSD map, possibly with no strain data  

%% change gID_list to 'ID_in_struCell'
gID_list = ID_in_struCell;
ind = ismember(gID,gID_list);
figure;
gd = gDiameter(ind);
max_gd = max(gd);
histogram(gd,0:10:120);
set(gca,'fontsize',18);
xlabel('Grain Diameter (\mum)');
ylabel('Counts');
print(fullfile(output_dir, 'Mg4Al grain diameter distribution'),'-dtiff');

%% Plot strain maps and for each iE, narrow colorbar range for exx, to differentiate tensile and compressive 
for iE = 1:6
    title_str{iE} = ['\fontsize{18}load step = ',num2str(iE),', \epsilon\fontsize{16}^G\fontsize{18} = ',num2str(strainPauses(iE),'%.4f')];
end

for iE = 1:6
    [f,a,c] = myplot(X,Y,dicData{iE}.exx,boundaryTFB);
    set(gca,'XTick',[],'YTick',[],'fontsize',18);
    title(c,'\fontsize{18}\epsilon_x_x');
    title(title_str{iE},'fontweight','normal');
    caxis([-0.01, 0.01]);
    print(fullfile(output_dir, ['color_adjusted_exx_',num2str(iE)]),'-dtiff');
    close;
end

%% special, strain map at iE = 3
iE = 4;
[f,a,c] = myplot(X,Y,dicData{iE}.exx,boundaryTFB);
caxis([-0.01, 0.01]);
set(gca,'XTick',[],'YTick',[],'fontsize',32, 'xlim',[37000, 47000], 'ylim',[12000,22000]);
set(c,'Ticks',[-0.01:0.005:0.1]);
print(fullfile(output_dir, ['color_adjusted_exx_',num2str(iE),'_zoom']),'-dtiff');
close;

%% plot ID map to help analysis
myplot(X,Y,ID,boundaryTFB);
label_map_with_ID(X,Y,ID,gcf,unique(ID(:)));


%% Manually select/change iE for analysis (iE = 1,2,3)
for iE = 1
    %%
    if iE == 1
        [f,a,c] = myplot(X, Y, dicData{iE}.exx, boundaryTFB);
    else
        [f,a,c] = myplot(X, Y, eEff{iE}, boundaryTFB);
    end
    label_map_with_trace(X,Y,ID, gID_list, 1, gca);     % label basal for all grains in AOI
    label_map_with_trace(X,Y,ID, 203, 5, gca);           % can label a special one with slip system other than basal
    title('double click to add grain, click x to exit', 'fontweight','normal');
    
    % first try to load previous result
    
    try
        load(fullfile(saveDataPath,['Mg4Al_gID_list_slip_iE_',num2str(iE),'.mat']), 'gID_list_slip', 'pos_list');
    catch
        if iE==1
            pos_list = [];  % [xpos1, ypos1; xpos2 ypos2; ...]
        else
            load([saveDataPath,'\Mg4Al_gID_list_slip_iE_',num2str(iE-1),'.mat'], 'gID_list_slip', 'pos_list', 'gID_list_special');
        end
    end
    
    
    gID_list_slip = []; % [grain_ID, slip_sys; ...]
    
    % draw existing ones, then keep select new ones, until hit 'x' to exit
    iPos = 1;
    while true
        try
            pos = pos_list(iPos,:);
            drawpoint(gca,'Position',pos);
        catch
            h = drawpoint;
            customWait(h);
            pos = round(h.Position);
        end
        [~,indc] = min(abs(pos(1)-X(1,:)));
        [~,indr] = min(abs(pos(2)-Y(:,1)));
        pos_list(iPos,:) = pos;
        gID_list_slip(iPos, :) = [ID(indr,indc), 1]     % by default, mark as basal slip.  Can modify later
        iPos = iPos + 1;
    end
    
    %% manually add [grain ID, slip sys] for non-basal;
    gID_list_special{1} = [203, 5];
    gID_list_special{2} = [203, 5;
        232, 5];
    gID_list_special{3} = [203, 5;
        232, 5];
    
    gID_list_slip = unique([gID_list_slip; gID_list_special{iE}], 'rows');
    save(fullfile(output_dir, ['Mg4Al_gID_list_slip_iE_',num2str(iE),'.mat']), 'gID_list_slip', 'pos_list', 'gID_list_special');       % can choose to save this for record
    
    %% [plot] label strain map with slip trace for selected grains
    % Regarding contrast, exx good for iE=1. eeff good for other iEs.
    if iE==1
        [f,a,c] = myplot(X,Y,dicData{iE}.exx,boundaryTFB);
        title(c,'\fontsize{18}\epsilon_x_x');
    else
        [f,a,c] = myplot(X,Y,eEff{iE},boundaryTFB);
        title(c,'\fontsize{18}\epsilon_e_f_f');
    end
    set(gca,'XTick',[],'YTick',[],'fontsize',18);
    title(title_str{iE},'fontweight','normal');
    
    % for each slip system, label grains
    unique_sss = unique(gID_list_slip(:,2));
    for ii = 1:length(unique_sss)
        ss_n = unique_sss(ii);
        inds = gID_list_slip(:,2)==ss_n;
        label_map_with_trace_for_Mg4Al(X,Y,ID, gID_list_slip(inds,1), ss_n, gca);
    end
    
    print(fullfile(output_dir, ['strain_map_with_trace_iE_',num2str(iE),'_with_trace']),'-dtiff');
    
    
end

%% Summarize slip and twin activity, grain wise and variant wise
% (1) summary for each grain
variableNames = {'ID', 'iE_initial_slip','max_basal_SF', 'iE_initial_twin','max_twin_SF'};
T1 = cell2table(cell(0,length(variableNames)));
T1.Properties.VariableNames = variableNames;

for ii = 1:length(gID_list)
    ID_current = gID_list(ii);

    ind = find(ID_current == gID);
    euler = [gPhi1(ind), gPhi(ind), gPhi2(ind)];
    [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler, [0 0 0], [0 0 0], stressTensor, sampleMaterial, 'twin');
    
    max_basal_SF = max(abs_schmid_factor(1:3,2));
                
    variant_SF = abs_schmid_factor(19:24,2);
    max_twin_SF = max(variant_SF);
    
    % initialize, but record max_basal_SF, max_twin_SF
    T1 = [T1; {ID_current, inf, max_basal_SF, inf, max_twin_SF}]; 
end

for iE = 1:3
    load(fullfile(output_dir,['Mg4Al_gID_list_slip_iE_',num2str(iE),'.mat']), 'gID_list_slip');      
    % check each ID, if has basal slip
    for ii = 1:length(gID_list)
        ID_current = gID_list(ii);        
        
        % If this grain has slip at iE
        ind = find(gID_list_slip(:,1) == ID_current);   % row_number or empty
        slip_sys = gID_list_slip(ind,2);    % active_slip_sys_numver or empty
        
        if ~isempty(slip_sys) && ismember(1,slip_sys) && isinf(T1.iE_initial_slip(ii))
            T1.iE_initial_slip(ii) = iE;
        end
        
        % If has twin active
        iS = find([struCell{iE}.gID] == ID_current);
        v_active = sum(struCell{iE}(iS).cTrueTwin,1)>0;
        
        if any(v_active) && isinf(T1.iE_initial_twin(ii))
            T1.iE_initial_twin(ii) = iE;
        end
    end
end

%  (2) summary for each variant
variableNames = {'ID', 'iE_initial_basal','max_basal_SF', 'iTwin','iE_initial_twin','variant_SF'};
T2 = cell2table(cell(0,length(variableNames)));
T2.Properties.VariableNames = variableNames;

for ii = 1:length(gID_list)
    ID_current = gID_list(ii);

    ind = find(ID_current == gID);
    euler = [gPhi1(ind), gPhi(ind), gPhi2(ind)];
    [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler, [0 0 0], [0 0 0], stressTensor, sampleMaterial, 'twin');
    
    max_basal_SF = max(abs_schmid_factor(1:3,2));
                
    variant_SF = abs_schmid_factor(19:24,2);
    max_twin_SF = max(variant_SF);
    
    ind = find(T1.ID==ID_current);
    iE_initial_slip = T1.iE_initial_slip(ind);
    
    % initialize, copy iE_initial_slip from T1, set iE_initial_twin as inf, record variant_SF, 
    for iTwin = 1:6
        T2 = [T2; {ID_current, iE_initial_slip, max_basal_SF, iTwin, inf, variant_SF(iTwin)}];
    end    
end

% for each iE, each gID, check each variant if active
for iE = 1:6
    twin_map = trueTwinMapCell{iE};
    for ii = 1:length(gID_list)
        ID_current = gID_list(ii);
        
        iS = find([struCell{iE}.gID] == ID_current);
        v_active = sum(struCell{iE}(iS).cTrueTwin,1)>0;
        
        % The above and the following are the same
        % inds = ismember(ID, ID_current);
        % variants = unique(twin_map(inds));
        % v_active = ismember(1:6, variants);
        
        for iTwin=1:6
           if v_active(iTwin)==1
              ind = (T2.ID==ID_current) & (T2.iTwin==iTwin);
              if isinf(T2.iE_initial_twin(ind))
                  T2.iE_initial_twin(ind) = iE;
              end
           end
        end        
    end
end


%% [plot scatterhist], show distributin of max_basal_SF vs max_twin_SF for all grains.
figure;
h = scatterhist(T1.max_twin_SF, T1.max_basal_SF, 'Marker', '.', 'MarkerSize', 14);
xlabel('Max Twin Schmid Factor');
ylabel('Max Basal Schmid Factor');
set(gca,'fontsize',14, 'xTick', [-0.1:0.1:0.5], 'yTick', [0:0.1:0.5]);
print(fullfile(output_dir, ['max basal SF vs max twin SF.tiff']), '-dtiff');

%% [plot eEff with trace] at iE=3, on eEff map, label slip systems. Label the ones shown slip in iE with '*' symbol === temp
iE = 3;
[f,a,c] = myplot(X, Y, eEff{iE}, boundaryTFB);
title(c,'\fontsize{18}\epsilon_e_f_f');
set(gca,'XTick',[],'YTick',[],'fontsize',18);
title(title_str{iE},'fontweight','normal');
caxis([0 0.09]);

% Grains with slip in iE=1
iE = 1;
load(fullfile(output_dir, ['Mg4Al_gID_list_slip_iE_',num2str(iE),'.mat']), 'gID_list_slip');
gID_list_slip_1 = gID_list_slip;

unique_sss = unique(gID_list_slip_1(:,2));
for ii = 1:length(unique_sss)
    ss_n = unique_sss(ii);
    inds = gID_list_slip_1(:,2)==ss_n;
    label_map_with_trace_for_Mg4Al(X,Y,ID, gID_list_slip_1(inds,1), ss_n, a);
end
label_map_with_text_a(X,Y,ID,f, 'target_ID', unique(gID_list_slip_1(:,1)), 'text', '1')

iE = 3;
load(fullfile(output_dir, ['Mg4Al_gID_list_slip_iE_',num2str(iE),'.mat']), 'gID_list_slip');
ind = ~ismember(gID_list_slip, gID_list_slip_1, 'rows');
gID_list_slip_3 = gID_list_slip(ind,:);

unique_sss = unique(gID_list_slip_3(:,2));
for ii = 1:length(unique_sss)
    ss_n = unique_sss(ii);
    inds = gID_list_slip_3(:,2)==ss_n;
    label_map_with_trace_for_Mg4Al(X,Y,ID, gID_list_slip_3(inds,1), ss_n, a);
end
label_map_with_text_a(X,Y,ID,f, 'target_ID', unique(gID_list_slip_3(:,1)), 'text', '3')

print(fullfile(output_dir, 'eeff map at iE=3 with traces.tiff'), '-dtiff');

%% [plot] pct of grains with observable basal slip traces at iE=1 and 3
edges = 0:0.05:0.5;
d_int = (edges(3)-edges(2))/2;
xpos = [edges(2)-d_int, edges(2:end-1)+d_int];
xstr = [];
for ii=1:length(edges)-1
    xstr{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
end
    
N_total = histcounts(T1.max_basal_SF, edges);

ind = T1.iE_initial_slip <= 1;    % iE=1
N_slip_1 = histcounts(T1.max_basal_SF(ind), edges);
pct_1 = N_slip_1./N_total;

ind = T1.iE_initial_slip <= 3;    % iE=3
N_slip_3 = histcounts(T1.max_basal_SF(ind), edges);
pct_3 = N_slip_3./N_total;

% plot
close;
figure; hold on;
plot(xpos, pct_1*100, '-o', 'color','r', 'linewidth',2);
plot(xpos, pct_3*100, '-s', 'color','b', 'linewidth',2);
set(gca,'fontsize',16,'xTick',xpos, 'xTickLabels',xstr,'xTickLabelRotation',45);
xlabel('Max Basal Schmid Factor');
ylabel('Percentage (%)');
legend('Load step 1', 'Load step 3', 'location','southeast');

print(fullfile(output_dir, 'pct grains slip iEs 1 3.tiff'), '-dtiff');













%% Analysis no longer needed:
%% Calculate grain Schmid factor. plot SF map and histogram
assert(eulerAligned==1,'check if euler is aligned to sample coordinate');
ss = define_SS_cart(sampleMaterial,'twin');

[gSF_basal, ~, ~, ~, gSF_etwin, ~, ~, ~, ~, ~, SF, SFs]...
    = calculate_SFs(gPhi1, gPhi, gPhi2, ss, 0, 0, 0, 0, 0, 0, stressTensor, sampleMaterial);

gSF_basal_map = assign_field_to_cell(ID, gID, gSF_basal, zeros(size(gID)));
gSF_etwin_map = assign_field_to_cell(ID, gID, gSF_etwin, zeros(size(gID)));

% [Figure] grain basal SF map
myplot(X,Y,gSF_basal_map,boundaryTFB);
set(gca,'XTick',[],'YTick',[],'fontsize',18);
title('Max Basal Schmid Factor','fontweight','normal','fontsize',18);
caxis([0, 0.5]);
print(fullfile(output_dir, ['basal_SF_map']),'-dtiff');

% [Figure] grain twin SF map
myplot(X,Y,gSF_etwin_map,boundaryTFB);
set(gca,'XTick',[],'YTick',[],'fontsize',18);
title('Max Twin Schmid Factor','fontweight','normal','fontsize',18);
caxis([-0.1, 0.5]);
print(fullfile(output_dir, ['etwin_SF_map']),'-dtiff');

% [Figure] grain basal SF histogram
figure;
t = gSF_basal(ismember(gID, gID_list));
histogram(t, 0:0.05:0.5);
ylabel('Counts');
xlabel('Max Basal Schmid Factor');
set(gca,'fontsize',18);
print(fullfile(output_dir, 'basal SF distribution'),'-dtiff');

% [Figure] grain twin SF histogram
figure;
t = gSF_etwin(ismember(gID, gID_list));
histogram(t, -0.1:0.05:0.5);
ylabel('Counts');
xlabel('Max Twin Schmid Factor');
set(gca,'fontsize',18);
print(fullfile(output_dir, 'twin SF distribution'),'-dtiff');

close all;



%% Schmid factor distribution of active slip system 
% ==> select iE = 1 or 2,3, and calculate the twin Schmid factor distribution of active twin 

iE = 3;
load([saveDataPath,'\Mg4Al_gID_list_slip_iE_3.mat'], 'gID_list_slip');  

% grain, max basal SF vs slip/non-slip
T1 = cell2table(cell(0,3));
T1.Properties.VariableNames = {'gID','SF_basal','slip_TF'};
t1_template = T1;

% variant, twin SF vs twin/non-twin
T2 = cell2table(cell(0,4));
T2.Properties.VariableNames = {'gID','variant','SF_twin','twin_TF'};
t2_template = T2;

% grain, max twin SF vs twin/non-twin
T3 = cell2table(cell(0,3));
T3.Properties.VariableNames = {'gID','SF_twin','twin_TF'};
t3_template = T3;


gID_list_twin = [];

warning off;
for iS = 1:length(struCell{iE})
    t1 = t1_template;
    t2 = t2_template;
    t3 = t3_template;
    
    ID_current = struCell{iE}(iS).gID;
    
    ind = find(ID_current == gID);
    euler = [gPhi1(ind), gPhi(ind), gPhi2(ind)];
    
    [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler, [0 0 0], [0 0 0], stressTensor, sampleMaterial, 'twin');
    
    basal_sf = max(abs_schmid_factor(1:3,2));
    slip_TF = ismember(ID_current, gID_list_slip);
    
    variant_sf = abs_schmid_factor(19:24,2);
    twin_TF = sum(struCell{iE}(iS).cTrueTwin,1)>0;
    
    if any(twin_TF)
       gID_list_twin = [gID_list_twin, ID_current];     % record twinned grain list  
    end
    
    t1.gID(1) = ID_current;
    t1.SF_basal(1) = basal_sf;
    t1.slip_TF(1) = slip_TF;
    T1 = [T1; t1];
    
    for ii = 1:6
        t2.gID(ii) = ID_current;
        t2.variant(ii) = ii;
        t2.SF_twin(ii) = variant_sf(ii);
        t2.twin_TF(ii) = twin_TF(ii);
    end
    T2 = [T2;t2];
    
    t3.gID(1) = ID_current;
    t3.SF_twin(1) = max(variant_sf);
    t3.twin_TF(1) = any(twin_TF);
    T3 = [T3; t3];
    
end
warning on;

%% double check, illustrate which grains twinned, which grains not twinned 
% map_notwin = trueTwinMapCell{iE};
% map_notwin(ismember(ID,gID_list_twin)) = -1; % make twinned grain -1
% myplot(map_notwin);
% 
% map_twin = trueTwinMapCell{iE};
% map_twin(~ismember(ID,gID_list_twin)) = -1;
% myplot(map_twin);

disp(['total # grains: ', num2str(length(gID_list))]);
disp(['total # grains twinned at iE = 3: ', num2str(length(gID_list_twin))]);
disp(['total # grains not twinned at iE = 3: ', num2str(length(gID_list)-length(gID_list_twin))]);

%% plot histograms to summarize distribution of SF, for slipped, not slipped, twinned, not twinned.  
edges_1 = 0:0.05:0.5;
edges_2 = [-0.3:0.05:0, 0.05:0.05:0.5];
edges_3 = [-0.1:0.05:0, 0.05:0.05:0.5];

% summarize basal_SF vs slip/non-slip
str = [];
for ii = 1:length(edges_1)-1
    str{ii} = [num2str(edges_1(ii)),'-',num2str(edges_1(ii+1))];
end
ind = T1.slip_TF==1;
h_T = histcounts(T1.SF_basal(ind), edges_1);
h_F = histcounts(T1.SF_basal(~ind), edges_1);
figure;
bar(edges_1(1:end-1), [h_T;h_F]');
set(gca,'XTick',edges_1(1:end-1), 'XTickLabel',str, 'XTickLabelRotation',45, 'fontsize',16);
xlabel('Max Basal Schmid Factor');
ylabel('Counts');
legend('Slip','No slip');
print('C:\Users\ZheChen\Desktop\slipTF vs basalSF iE select','-dtiff');

% summarize twin_SF vs twin/no-twin, variant-wise
str = [];
for ii = 1:length(edges_2)-1
    str{ii} = [num2str(edges_2(ii),2),'-',num2str(edges_2(ii+1),2)];
end
ind = T2.twin_TF==1;
h_T = histcounts(T2.SF_twin(ind), edges_2);
h_F = histcounts(T2.SF_twin(~ind), edges_2);
figure;
bar(edges_2(1:end-1), [h_T;h_F]');
set(gca,'XTick',edges_2(1:end-1), 'XTickLabel',str, 'XTickLabelRotation',45, 'fontsize',16);
xlabel('Twin Schmid Factor');
ylabel('Counts');
legend('Twinned','Not twinned');
print('C:\Users\ZheChen\Desktop\twinTF vs variantSF','-dtiff');

% summarize twin_SF vs slip/non-slip, grain-wise
str = [];
for ii = 1:length(edges_3)-1
    str{ii} = [num2str(edges_3(ii)),'-',num2str(edges_3(ii+1))];
end
ind = T3.twin_TF==1;
h_T = histcounts(T3.SF_twin(ind), edges_3);
h_F = histcounts(T3.SF_twin(~ind), edges_3);
figure;
bar(edges_3(1:end-1), [h_T;h_F]');
set(gca,'XTick',edges_3(1:end-1), 'XTickLabel',str, 'XTickLabelRotation',45, 'fontsize',16);
xlabel('Max Twin Schmid Factor');
ylabel('Counts');
legend('Twinned','Not twinned','location','northwest');
print('C:\Users\ZheChen\Desktop\twinTF vs twinSF','-dtiff');



























