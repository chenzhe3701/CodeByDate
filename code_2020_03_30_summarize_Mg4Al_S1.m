% based on script_PRISMS_2019_twinPct_with_error_bar
% chenzhe 2020-03-30
% Plot the evolution of twin area fraction vs global strain level.
% Divide the manually corrected twin map into 3x3=9 regions for errorbar

clear twinPct;

load('D:\p\m\DIC_Analysis\20200324_2004_relabeled_result_Mg4Al_S1.mat','trueTwinMapCell');
load('D:\p\m\DIC_Analysis\setting_for_real_samples\Mg4Al_S1_setting.mat','strainPauses', 'strainPauses_sg');
load('E:\Mg4Al_S1_insitu\Analysis_by_Matlab\Mg4Al_S1_EbsdToSemForTraceAnalysis.mat','ID','boundaryTFB');
load('D:\p\m\DIC_Analysis\temp_results\20200325_0506_new_variant_map_Mg4Al_S1.mat','variantMapCleanedCell');
strainFileFolder = 'E:\Mg4Al_S1_insitu\SEM_Images\stitched_DIC';
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



