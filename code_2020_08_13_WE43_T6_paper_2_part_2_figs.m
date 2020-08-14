
clear;
clc;
% settings and paths
[fileSetting,pathSetting] = uigetfile('D:\p\m\DIC_Analysis\setting_for_real_samples\WE43_T6_C1_setting.mat','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','sampleMaterial','stressTensor','strainPauses');

saveDataPath = 'D:\WE43_T6_C1\Analysis_by_Matlab_after_realign';
dicPath = 'D:\WE43_T6_C1\SEM Data\stitched_DIC';
variantFile = 'D:\p\m\DIC_Analysis\temp_results\WE43_T6_C1_new_variant_map_20200401.mat';

for iE = 1:5
    dicFileName = ['_',num2str(iE),'_v73.mat'];
    dicData{iE} = matfile(fullfile(dicPath, dicFileName));
end

load(fullfile(saveDataPath, 'WE43_T6_C1_EbsdToSemForTraceAnalysis_GbAdjusted.mat'),...
    'ID','X','Y','boundaryTF','boundaryTFB','gDiameter','gID','gNeighbors','gNNeighbors','gPhi1','gPhi','gPhi2','eulerAligned');
load(variantFile,'variantMapCleanedCell','struCell');
%%
for iE = 1:5
    strain_str{iE} = ['\epsilon\fontsize{16}^G = ',num2str(strainPauses(iE),'%.3f')];
end

clims = [-0.014    0.014
   -0.03    0.014
   -0.06    0.014
   -0.08    0.014
   -0.10    0.014];  % climits for each strain level
for iE = 1:5
    %     eEff{iE} = calculate_effective_strain(dicData{iE}.exx, dicData{iE}.exy, dicData{iE}.eyy);
    %     myplot(X,Y,eEff{iE},boundaryTFB);
    exx = dicData{iE}.exx;
    exx((exx>0.017)&(X>28500)&(X<32500)&(Y<3500)) = nan;
    [f,a,c] = myplot(X,Y,exx,boundaryTFB);
    %     drawline('position',[1000,35000; 1000 + 2276*5, 35000], 'linewidth',100);
    set(gca,'xTick',[],'yTick',[],'fontsize',18);
    if iE==1
       set(c,'TickLabels',{'-0.01','','0','','0.01'});
    end
    caxis(clims(iE,:));
%     clims(iE,:) = caxis;
    title(strain_str(iE),'fontweight','normal');
    title(c,'\fontsize{20}\epsilon\fontsize{18}_x_x');
    print(['C:\Users\ZheChen\Desktop\exxMap_',num2str(iE)],'-r300','-dtiff');
end

% for iE = 3
% 
%     myplot(X,Y,dicData{iE}.exx > 0.0,boundaryTFB);
% 
% end

for iE = 2:5
    [f,a,c] = myplot(X,Y,variantMapCleanedCell{iE},boundaryTFB);
    caxis([-0.5,6.5]);
    cmap = parula(6);
    cmap = jet(6);
    cmap = [0.3, 0.3, 0.3; cmap];
    colormap(cmap);
    set(c,'limits', [0.5, 6.5]);
    set(gca,'XTick',[],'YTick',[],'fontsize',18);
    title(strain_str(iE),'fontweight','normal');
    print(['C:\Users\ZheChen\Desktop\variantMap_',num2str(iE)],'-r300','-dtiff');
end

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
print(['C:\Users\ZheChen\Desktop\basal_SF_map'],'-r300','-dtiff');

% [Figure] grain twin SF map
myplot(X,Y,gSF_etwin_map,boundaryTFB);
set(gca,'XTick',[],'YTick',[],'fontsize',18);
title('Max Twin Schmid Factor','fontweight','normal','fontsize',18);
caxis([-0.1, 0.5]);
print(['C:\Users\ZheChen\Desktop\etwin_SF_map'],'-r300','-dtiff');

gID_list = unique(ID(:));
% [Figure] grain basal SF histogram
figure;
t = gSF_basal(ismember(gID, gID_list));
histogram(t, 0:0.05:0.5);
ylabel('Counts');
xlabel('Max Basal Schmid Factor');
set(gca,'fontsize',18);
print('C:\Users\ZheChen\Desktop\basal SF distribution','-dtiff');

% [Figure] grain twin SF histogram
figure;
t = gSF_etwin(ismember(gID, gID_list));
histogram(t, -0.1:0.05:0.5);
ylabel('Counts');
xlabel('Max Twin Schmid Factor');
set(gca,'fontsize',18);
print('C:\Users\ZheChen\Desktop\twin SF distribution','-dtiff');

%%

