
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

% twin variant data
variant_file = 'E:\Mg4Al_S1_insitu\Analysis_by_Matlab\previous result recovered\20200325_0506_new_variant_map_Mg4Al_S1.mat';
load(variant_file, 'struCell','variantMapCell');

output_dir = 'E:\zhec umich Drive\0_temp_output\Mg4Al_S1 analysis for paper b';
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

ind = ismember(gID,gID_list);
figure;
gd = gDiameter(ind);
max_gd = max(gd);
histogram(gd,0:10:120);
set(gca,'fontsize',18);
xlabel('Grain Diameter (\mum)');
ylabel('Counts');
print(fullfile(output_dir, '\Mg4Al grain diameter distribution'),'-dtiff');

%% Plot strain maps
for iE = 1:6
   strain_str{iE} = ['\epsilon\fontsize{16}^G = ',num2str(strainPauses(iE),'%.4f')]; 
end
for iE = 1:6
    [f,a,c] = myplot(X,Y,dicData{iE}.exx,boundaryTFB);
    set(gca,'XTick',[],'YTick',[],'fontsize',18);
    caxis([-0.08, 0.01]);
    title(c,'\fontsize{18}\epsilon_x_x');
    title(strain_str{iE},'fontweight','normal');
    print(fullfile(output_dir, ['exx_',num2str(iE)]),'-dtiff');
end
close all;

%% For each iE, narrow colorbar range to differentiate tensile and compressive 
for iE = 1:6
    [f,a,c] = myplot(X,Y,dicData{iE}.exx,boundaryTFB);
    set(gca,'XTick',[],'YTick',[],'fontsize',18);
    title(c,'\fontsize{18}\epsilon_x_x');
    title(strain_str{iE},'fontweight','normal');
    caxis([-0.01, 0.01]);
    print(fullfile(output_dir,['color_adjusted_exx_',num2str(iE)]),'-dtiff');
end
close all;

%% Plot twin variant maps
for iE = 1:6
   eG_str{iE} = ['\epsilon\fontsize{16}^G = ',num2str(strainPauses(iE),'%.4f')]; 
end
for iE = 1:6
    [f,a,c] = myplot(X,Y,variantMapCell{iE},boundaryTFB);
    caxis([-0.5,6.5]);
    cmap = parula(6);
    cmap = [0.3, 0.3, 0.3; cmap];
    colormap(cmap);
    set(c,'limits', [0.5, 6.5]);
    set(gca,'XTick',[],'YTick',[],'fontsize',18);
    title(eG_str{iE},'fontweight','normal');
    print(fullfile(output_dir,['variantMap_',num2str(iE)]),'-dtiff');
end
close all;

%% plot ID map
myplot(X,Y,ID,boundaryTFB);
label_map_with_ID(X,Y,ID,gcf,unique(ID(:)));

%% (step 1) @ iE=1, use exx_map to actively selected grains with slip trace.  
addpath('D:\p\m\DIC_Analysis\on_going_work');

iE = 1;
[f,a,c] = myplot(X,Y,dicData{iE}.exx,boundaryTFB);
label_map_with_trace(X,Y,ID,gID_list, 1, gca);
label_map_with_trace(X,Y,ID,203, 5, gca);
title('double click to add grain, click x to exit', 'fontweight','normal'); 

try
    load([saveDataPath,'\Mg4Al_gID_list_slip_iE_1.mat'], 'gID_list_slip', 'pos_list', 'gID_list_special');       % first try to load previous result
catch
    pos_list = [];
end
gID_list_slip = [];    

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
    gID_list_slip(iPos) = ID(indr,indc)
    iPos = iPos + 1;
end
%%
gID_list_special = [203, 5];    % manually add [grain ID, slip system other than basal];

if ~exist(fullfile(saveDataPath,'Mg4Al_gID_list_slip_iE_1.mat'),'file')
    save(fullfile(saveDataPath,'Mg4Al_gID_list_slip_iE_1.mat'), 'gID_list_slip', 'pos_list', 'gID_list_special');       % can choose to save this for record
end
%% (step 2) label strain map with slip trace for selected grains at iE=1
iE = 1;
[f,a,c] = myplot(X,Y,dicData{iE}.exx,boundaryTFB);
set(gca,'XTick',[],'YTick',[],'fontsize',18);

label_map_with_trace_for_Mg4Al(X,Y,ID,gID_list_slip, 1, gca);
for ii = 1:size(gID_list_special,1)
    label_map_with_trace_for_Mg4Al(X,Y,ID, gID_list_special(ii,1), gID_list_special(ii,2), gca);   % can modify this code to plot trace as black
end

title(c,'\fontsize{18}\epsilon_x_x');
title(strain_str{iE},'fontweight','normal');
print(fullfile(output_dir, ['exx_',num2str(iE),'_with_trace']),'-dtiff');
close all;

%% (step 3) at iE = 2 (or higher, can change), continue actively selected grains with slip trace, with eEff map  
iE = 2;
[f,a,c] = myplot(X,Y,eEff{iE},boundaryTFB);

label_map_with_trace(X,Y,ID,gID_list, 1, gca);
label_map_with_trace(X,Y,ID,203, 5, gca);   % something special to be modified   

try
    load([saveDataPath,'\Mg4Al_gID_list_slip_iE_2.mat'], 'gID_list_slip', 'pos_list', 'gID_list_special');       % first try to load previous result
catch
    load([saveDataPath,'\Mg4Al_gID_list_slip_iE_1.mat'], 'gID_list_slip', 'pos_list', 'gID_list_special');       % first try to load previous result
end
gID_list_slip = [];    

% draw existing ones, then keep select new ones, until hit 'x' to exit
iPos = 1;
while true
    try
        pos = pos_list(iPos,:);
        drawpoint(gca,'Position',pos,'color','r');
    catch
        h = drawpoint;
        customWait(h);
        pos = round(h.Position);
    end
    [~,indc] = min(abs(pos(1)-X(1,:)));
    [~,indr] = min(abs(pos(2)-Y(:,1)));
    pos_list(iPos,:) = pos;
    gID_list_slip(iPos) = ID(indr,indc)
    iPos = iPos + 1;
end
%%
gID_list_special = [203, 5;
    232, 5];    % manually add [grain ID, slip system other than basal]; -> grain 52, don't know ss.
if ~exist(fullfile(saveDataPath,'Mg4Al_gID_list_slip_iE_2.mat'),'file')
    save(fullfile(saveDataPath,'Mg4Al_gID_list_slip_iE_2.mat'), 'gID_list_slip', 'pos_list', 'gID_list_special');       % can choose to save this for record
end
%% (step 4) label strain map with slip trace for selected grains at iE=2
iE = 2;
[f,a,c] = myplot(X,Y,eEff{iE},boundaryTFB);
set(gca,'XTick',[],'YTick',[],'fontsize',18);

label_map_with_trace_for_Mg4Al(X,Y,ID,gID_list_slip, 1, gca);
for ii = 1:size(gID_list_special,1)
    label_map_with_trace_for_Mg4Al(X,Y,ID, gID_list_special(ii,1), gID_list_special(ii,2), gca);   % can modify this code to plot trace as black
end

title(c,'\fontsize{18}\epsilon_e_f_f');
title(strain_str{iE},'fontweight','normal');
print(fullfile(output_dir, ['eeff_',num2str(iE),'_with_trace']),'-dtiff');
close all;

%% (step 5) at iE = 3 (or higher, can change), continue actively selected grains with slip trace, with eEff map  
iE = 3;
[f,a,c] = myplot(X,Y,eEff{iE},boundaryTFB);

label_map_with_trace(X,Y,ID,gID_list, 1, gca);
label_map_with_trace(X,Y,ID,203, 5, gca);   % something special to be modified   

try
    load([saveDataPath,'\Mg4Al_gID_list_slip_iE_3.mat'], 'gID_list_slip', 'pos_list', 'gID_list_special');       % first try to load previous result
catch
    load([saveDataPath,'\Mg4Al_gID_list_slip_iE_2.mat'], 'gID_list_slip', 'pos_list', 'gID_list_special');       % first try to load previous result
end
gID_list_slip = [];    

% draw existing ones, then keep select new ones, until hit 'x' to exit
iPos = 1;
while true
    try
        pos = pos_list(iPos,:);
        drawpoint(gca,'Position',pos,'color','r');
    catch
        h = drawpoint;
        customWait(h);
        pos = round(h.Position);
    end
    [~,indc] = min(abs(pos(1)-X(1,:)));
    [~,indr] = min(abs(pos(2)-Y(:,1)));
    pos_list(iPos,:) = pos;
    gID_list_slip(iPos) = ID(indr,indc)
    iPos = iPos + 1;
end
%%
gID_list_special = [203, 5;
    232, 5];    % manually add [grain ID, slip system other than basal]; -> grain 52, don't know ss.

if ~exist(fullfile(saveDataPath,'Mg4Al_gID_list_slip_iE_3.mat'),'file')
    save(fullfile(saveDataPath,'Mg4Al_gID_list_slip_iE_3.mat'), 'gID_list_slip', 'pos_list', 'gID_list_special');       % can choose to save this for record
end
%% (step 6) label strain map with slip trace for selected grains at iE=2
iE = 3;
[f,a,c] = myplot(X,Y,eEff{iE},boundaryTFB);
set(gca,'XTick',[],'YTick',[],'fontsize',18);

label_map_with_trace_for_Mg4Al(X,Y,ID,gID_list_slip(~ismember(gID_list_slip, gID_list_special(:,1))), 1, gca);
for ii = 1:size(gID_list_special,1)
    label_map_with_trace_for_Mg4Al(X,Y,ID, gID_list_special(ii,1), gID_list_special(ii,2), gca);   % can modify this code to plot trace as black
end

title(c,'\fontsize{18}\epsilon_e_f_f');
title(strain_str{iE},'fontweight','normal');
print(fullfile(output_dir,['eeff_',num2str(iE),'_with_trace']),'-dtiff');

close all;

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
print(fullfile(output_dir, 'basal_SF_map'),'-dtiff');

% [Figure] grain twin SF map
myplot(X,Y,gSF_etwin_map,boundaryTFB);
set(gca,'XTick',[],'YTick',[],'fontsize',18);
title('Max Twin Schmid Factor','fontweight','normal','fontsize',18);
caxis([-0.1, 0.5]);
print(fullfile(output_dir, 'etwin_SF_map'),'-dtiff');

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
print(fullfile(output_dir,'twin SF distribution'),'-dtiff');

close all;

%% [] explore. struCell has only 177 grains, gID_list has 180. find the difference
ID_in_struCell = [];
for ii = 1:length(struCell{2})
    ID_in_struCell = [ID_in_struCell, struCell{2}(ii).gID];
end
ID_diff = gID_list(~ismember(gID_list, ID_in_struCell));
myplot(ismember(ID,ID_diff));
gDiameter(ismember(gID,ID_diff))
% -> Note: the result shows that the 3 grains are edge grains on EBSD map, possibly with no strain data  

%% Calculate the slip Schmid factor distribution of active slip system at iE = 1 or 2,3 --> select
%% And calculate the twin Schmid factor distribution of active twin 

for iE = 1:3
load([saveDataPath,'\Mg4Al_gID_list_slip_iE_',num2str(iE),'.mat'], 'gID_list_slip');  

T1 = cell2table(cell(0,3));
T1.Properties.VariableNames = {'gID','SF_basal','slip_TF'};
t1_template = T1;

T2 = cell2table(cell(0,4));
T2.Properties.VariableNames = {'gID','variant','SF_twin','twin_TF'};
t2_template = T2;

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
map_notwin = variantMapCell{iE};
map_notwin(ismember(ID,gID_list_twin)) = -1; % make twinned grain -1
myplot(map_notwin);

map_twin = variantMapCell{iE};
map_twin(~ismember(ID,gID_list_twin)) = -1;
myplot(map_twin);

disp(['total # grains: ', num2str(length(gID_list))]);
disp(['total # grains twinned at iE = 3: ', num2str(length(gID_list_twin))]);
disp(['total # grains not twinned at iE = 3: ', num2str(length(gID_list)-length(gID_list_twin))]);

%% plot histograms to summarize distribution of SF, for slipped, not slipped, twinned, not twinned.  
close all;
edges_basal = 0:0.05:0.5;
edges_2 = [-0.3:0.05:0, 0.05:0.05:0.5];
edges_3 = [-0.1:0.05:0, 0.05:0.05:0.5];

% summarize basal_SF vs slip/non-slip
str = [];
for ii = 1:length(edges_basal)-1
    str{ii} = [num2str(edges_basal(ii)),'-',num2str(edges_basal(ii+1))];
end
ind1 = T1.slip_TF==1;
ind2 = T1.slip_TF==0;
h_T = histcounts(T1.SF_basal(ind1), edges_basal);    % with slip trace
h_F = histcounts(T1.SF_basal(ind2), edges_basal);   % without slip trace
figure;
hbar = bar(edges_basal(1:end-1)+0.025, [h_F;h_T]', 0.9, 'stacked');
set(gca,'XTick',edges_basal(1:end-1)+0.025, 'XTickLabel',str, 'XTickLabelRotation',45, 'fontsize',16);
xlabel('Max Basal Schmid Factor');
ylabel('Counts');
hbar(1).FaceColor = [0 0 1];
hbar(2).FaceColor = [1 0 0];
    
yyaxis right;
plot(edges_basal(1:end-1)+0.025, h_T./(h_T+h_F)*100, '-ok', 'linewidth',1.5);   
set(gca,'ycolor','k','ylim',[0 160]);    
ylabel('Percent (%)');

legend('Without Slip Traces','With Slip Traces', 'Percent with Slip Traces');
print(fullfile(output_dir, ['slipTF vs basalSF iE=',num2str(iE),'.tiff']),'-dtiff');


% summarize twin_SF vs twin/no-twin, variant-wise
str = [];
for ii = 1:length(edges_2)-1
    str{ii} = [num2str(edges_2(ii),2),'-',num2str(edges_2(ii+1),2)];
end
ind1 = T2.twin_TF==1;
ind2 = T2.twin_TF==0;
h_T = histcounts(T2.SF_twin(ind1), edges_2);
h_F = histcounts(T2.SF_twin(ind2), edges_2);
figure;
bar(edges_2(1:end-1), [h_T;h_F]');
set(gca,'XTick',edges_2(1:end-1), 'XTickLabel',str, 'XTickLabelRotation',45, 'fontsize',16);
xlabel('Twin Schmid Factor');
ylabel('Counts');
legend('Twinned','Not twinned');
% print('C:\Users\ZheChen\Desktop\twinTF vs variantSF','-dtiff');


% summarize twin_SF vs slip/non-slip, grain-wise
str = [];
for ii = 1:length(edges_3)-1
    str{ii} = [num2str(edges_3(ii)),'-',num2str(edges_3(ii+1))];
end
ind1 = T3.twin_TF==1;
ind2 = T3.twin_TF==0;
h_T = histcounts(T3.SF_twin(ind1), edges_3);
h_F = histcounts(T3.SF_twin(ind2), edges_3);
figure;
bar(edges_3(1:end-1), [h_T;h_F]');

set(gca,'XTick',edges_3(1:end-1), 'XTickLabel',str, 'XTickLabelRotation',45, 'fontsize',16);
xlabel('Max Twin Schmid Factor');
ylabel('Counts');
legend('Twinned','Not twinned','location','northwest');
% print('C:\Users\ZheChen\Desktop\twinTF vs twinSF','-dtiff');

%% double check, illustrate which grains do not have slip traces
map_no_slip = ID;
map_no_slip(ismember(ID, gID_list_slip)) = 0;
myplot(map_no_slip,boundaryTFB);
print(fullfile(output_dir, ['grains no basal slip trace iE=',num2str(iE),'.tiff']),'-dtiff');
end

























