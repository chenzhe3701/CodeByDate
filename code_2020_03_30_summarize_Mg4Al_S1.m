% based on script_PRISMS_2019_twinPct_with_error_bar
% chenzhe 2020-03-30
% Plot the evolution of twin area fraction vs global strain level.
% Divide the manually corrected twin map into 3x3=9 regions for errorbar

clear twinPct;
for_TMS_dir = 'E:\zhec umich Drive\0_temp_output';
save_dir = 'E:\Mg4Al_S1_insitu\Summary';

load('D:\p\m\DIC_Analysis\20200324_2004_relabeled_result_Mg4Al_S1.mat','trueTwinMapCell');
load('D:\p\m\DIC_Analysis\setting_for_real_samples\Mg4Al_S1_setting.mat','strainPauses', 'strainPauses_sg');
load('E:\Mg4Al_S1_insitu\Analysis_by_Matlab\Mg4Al_S1_EbsdToSemForTraceAnalysis.mat','ID','boundaryTFB');
load('D:\p\m\DIC_Analysis\temp_results\20200325_0506_new_variant_map_Mg4Al_S1.mat','variantMapCleanedCell');
strainFileFolder = 'E:\Mg4Al_S1_insitu\SEM Data\stitched_DIC';
iE_start = 1;
iE_stop = 6;
for iE = iE_start:iE_stop
    variantMap = variantMapCleanedCell{iE};
    if iE==iE_start
        [nR,nC] = size(variantMap);
    end
    nr = 3;
    nc = 3;
    for ir=1:nr
       for ic = 1:nc
           subTrueTwinMap = variantMap([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic);
           subIDMap = ID([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic);
           twinPct((ir-1)*nc+ic,iE) = sum(subTrueTwinMap(:)>0)/sum(subIDMap(:)>0);
       end
    end
    
end

%%
close all;
tAvg = mean(twinPct);
tStd = std(twinPct);

% (1) DIC strain
figure;
eg(iE_start:iE_stop) = strainPauses(iE_start:iE_stop);    % eg(1:6) = [-0.0012, -0.0117, -0.0186, -0.00172, -0.0124, 0]; % this is from DIC strain

errorbar(eg(iE_start:iE_stop), 100*tAvg(iE_start:iE_stop), 100*tStd(iE_start:iE_stop), '-r.','linewidth',2,'markersize',18);
for ii=iE_start:iE_stop
    if ii==3
        text(strainPauses(ii)+0, 3+tAvg(ii)*100,[num2str(100*tAvg(ii),2),'\pm',num2str(100*tStd(ii),2),'%'],'fontsize',14);
    else
        text(strainPauses(ii)+0.001, 0+tAvg(ii)*100,[num2str(100*tAvg(ii),2),'\pm',num2str(100*tStd(ii),2),'%'],'fontsize',14);
    end
end
set(gca,'linewidth',1.5);
set(gca,'xlim',[-0.025, 0.01],'ylim',[-2 28],'fontsize',18,'fontweight','normal');
xlabel('\fontsize{24}\epsilon\fontsize{18}^G, DIC Strain');
ylabel('Twin Area Percent, %');
print('E:\Mg4Al_S1_insitu\Summary\Twin Area vs DIC Strain.tif','-dtiff');
% (2) strain gage strain
figure;
eg(iE_start:iE_stop) = strainPauses_sg(iE_start:iE_stop);    % eg(1:6) = [-0.0012, -0.0117, -0.0186, -0.00172, -0.0124, 0]; % this is from DIC strain

errorbar(eg(iE_start:iE_stop), 100*tAvg(iE_start:iE_stop), 100*tStd(iE_start:iE_stop), '-r.','linewidth',2,'markersize',18);
for ii=iE_start:iE_stop
    if ii==3
        text(eg(ii)+0, 3+tAvg(ii)*100,[num2str(100*tAvg(ii),2),'\pm',num2str(100*tStd(ii),2),'%'],'fontsize',14);
    else
        text(eg(ii)+0.001, 0+tAvg(ii)*100,[num2str(100*tAvg(ii),2),'\pm',num2str(100*tStd(ii),2),'%'],'fontsize',14);
    end
end
set(gca,'linewidth',1.5);
set(gca,'xlim',[-0.03, 0.015],'ylim',[-2 28],'fontsize',18,'fontweight','normal');
xlabel('\fontsize{24}\epsilon\fontsize{18}^G, Strain Gage Strain');
ylabel('Twin Area Percent, %');
print('E:\Mg4Al_S1_insitu\Summary\Twin Area vs Strain Gage Strain.tif','-dtiff');

%% Include iE=0 as tAvg{1}. 
tAvg(2:7) = tAvg(end-5:end);
tStd(2:7) = tStd(end-5:end);
twinPct(2:7) = twinPct(end-5:end);
save(fullfile(save_dir, 'twin_pct.mat'), 'twinPct', 'tAvg', 'tStd');

strain_dic = [0, -0.0012, -0.0117, -0.0186, ...
    -0.0172, -0.0124, 0.0006];

strain_sg = [0, -0.0075, -0.0150, -0.0250, ...
    -0.0230, -0.0170, 0];

colors = parula(5);

% [1] using strain gage strain
figure; hold on;
inds = {1:4, 4:7};

errorbar(strain_sg(inds{1}), 100*tAvg(inds{1}), 100*tStd(inds{1}), '.-', 'color',colors(1,:), 'linewidth',1.5,'markersize',24);
errorbar(strain_sg(inds{2}), 100*tAvg(inds{2}), 100*tStd(inds{2}), '.-', 'color',colors(2,:), 'linewidth',1.5,'markersize',24);

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',[-0.035, 0.005],'ylim',[-2 30],'fontsize',18,'fontweight','normal');
xlabel('Strain from strain gage');
ylabel('Twin Area Percent (%)');
print(fullfile(save_dir,'twin_pct_vs_sg.tiff'),'-dtiff');

% [2] using (fine transformed) ebsd estimated strain
figure; hold on;
inds = {1:4, 4:7};

errorbar(strain_dic(inds{1}), 100*tAvg(inds{1}), 100*tStd(inds{1}), '.-', 'color',colors(1,:), 'linewidth',1.5,'markersize',24);
errorbar(strain_dic(inds{2}), 100*tAvg(inds{2}), 100*tStd(inds{2}), '.-', 'color',colors(2,:), 'linewidth',1.5,'markersize',24);

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',[-0.025, 0.005],'ylim',[-2 30],'fontsize',18,'fontweight','normal');
xlabel('Strain from ebsd estimate');
ylabel('Twin Area Percent (%)');
print(fullfile(save_dir,'twin_pct_vs_dic_strain.tiff'),'-dtiff');


tbl = array2table([(0:6)', strain_sg(:), strain_dic(:), 100*tAvg(:), 100*tStd(:)]);
tbl.Properties.VariableNames = {'iE','strain_sg','strain_dic','twinPct %','twinStd %'};
disp(tbl);
save(fullfile(save_dir, 'twin_pct.mat'), 'twinPct', 'tAvg', 'tStd', 'tbl');
%%

% print('E:\Mg4Al_S1_insitu\Summary\Twin Area vs Global Strain.tif','-dtiff');

%% This is a version that gives a little bit different result, as it only calculates one average of the whole map
for iE=iE_start:iE_stop
    aF(iE) = sum(trueTwinMapCell{iE}(:)>0) / sum(ID(:)>0)
end

figure; 
plot(strainPauses(iE_start:iE_stop),aF(iE_start:iE_stop)*100,'-o');
set(gca,'fontsize',16)
xlabel('strain')
ylabel('twin area percentage')

for ii=iE_start:iE_stop
    if ii==3
        text(strainPauses(ii)+0, 3+tAvg(ii)*100,[num2str(100*aF(ii),2),'%'],'fontsize',14);
    else
        text(strainPauses(ii)+0.001, 0+tAvg(ii)*100,[num2str(100*aF(ii),2),'%'],'fontsize',14);
    end
end

%% Plot and save variantMap
for iE=1:6
   myplot(variantMapCleanedCell{iE},boundaryTFB);
   caxis([0 6]);
   print(['E:\Mg4Al_S1_insitu\Summary\variantMap_',num2str(iE),'.tiff'],'-dtiff');
end

%% plot the exx maps
for iE=1:6
    load(fullfile(strainFileFolder,['_',num2str(iE),'_v73.mat']),'exx');
    myplot(exx,boundaryTFB);
%     caxis([-0.1 0.02]);
   print(['E:\Mg4Al_S1_insitu\Summary\exxMap_',num2str(iE),'.tiff'],'-dtiff');
end

%% code
clear eulerAligned;
load('E:\Mg4Al_S1_insitu\Analysis_by_Matlab\Mg4Al_S1_EbsdToSemForTraceAnalysis.mat','ID','gID','gPhi1','gPhi','gPhi2','phi1','phi','phi2','eulerAligned');
if ~exist('eulerAligned','var')
    error('check eulerAligned before continue')
end

save('E:\Mg4Al_S1_insitu\Summary\Mg4Al_S1_EBSD_organized.mat','ID','gID','gPhi1','gPhi','gPhi2','phi1','phi','phi2','eulerAligned','-v7.3');

%% maybe can also choose to save variantMapCleanedCell
save('E:\Mg4Al_S1_insitu\Summary\twinVariantMapCell.mat','variantMapCleanedCell','-v7.3');

%% For TMS 2011, Plot and save variantMap, gray background
cmap = parula(6);
cmap = [[1,1,1] * 0.3; cmap];
for iE=1:6
   [f,a,c] = myplot(variantMapCleanedCell{iE},boundaryTFB);
   caxis([-0.5 6.5]);
   colormap(cmap);
   set(c,'limits',[0.5, 6.5]);
   axis off;
   print(fullfile(for_TMS_dir, ['Mg4Al_S1_variantMap_',num2str(iE),'.tiff']),'-dtiff');
   close;
end

%% For TMS 2011, plot the exx maps
for iE=1:6
    load(fullfile(strainFileFolder,['_',num2str(iE),'_v73.mat']),'exx','x','y');
    myplot(exx,boundaryTFB);
%     caxis([-0.1 0.02]);
    title('');
    axis off;
   print(fullfile(for_TMS_dir, ['Mg4Al_S1_exxMap_',num2str(iE),'.tiff']),'-dtiff');
   close;
end

