% Get grain data for grian 639, which is the grain used as example in 2020
% Mater Char paper.  make plots for TMS-2021 slides

clear;
addChenFunction;

% dir to save data
working_dir = 'E:\zhec umich Drive\0_temp_output\For 2021 TMS';

% DIC data path
dicPath = 'D:\WE43_T6_C1\SEM Data\stitched_DIC';
% sampleName = 'WE43_T6_C1';

% data for X, Y, ID, etc
saveDataPath = [uigetdir('D:\WE43_T6_C1\Analysis_by_Matlab_after_realign','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;

if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end
try
    load([saveDataPath,'WE43_T6_C1_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','cityDistMap','ID','gID','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');
catch
    load([saveDataPath,'WE43_T6_C1_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','cityDistMap','ID','gID','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');
end

% Identfied variant data
[newVariantFile, newVariantPath] = uigetfile('D:\p\m\DIC_Analysis\temp_results\WE43_T6_C1_new_variant_map.mat','select the results for variantMapCleanedCell');
load(fullfile(newVariantPath, newVariantFile),'variantMapCleanedCell','struCell');

% Make unique grain boundary map, and a list of the unique grain boundaries
[~, boundaryID, neighborID, ~, ~] = find_one_boundary_from_ID_matrix(ID);
uniqueBoundary = max(boundaryID,neighborID)*10000 + min(boundaryID,neighborID);
uniqueBoundaryList = unique(uniqueBoundary(:));
uniqueBoundaryList(uniqueBoundaryList==0)=[];

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

[ssa, c_a, nss, ntwin, ssGroup] = define_SS('Mg','twin');

% Make sure data is exy_corrected, and outlier_removed.  Convert into v7.3 for partial loading
clear strainFile;
for iE = iE_start:iE_stop
    strainFileName = [dicPath,'\','_',num2str(iE)];
    disp(strainFileName);
    if ~exist([strainFileName,'_v73.mat'],'file')
        load(strainFileName);
        clear('exy_corrected');
        load(strainFileName,'exy_corrected');   % if 'exy_corrected' does not exist, this does not give error, rather, just warning.
        
        if exist('exy_corrected','var')&&(1==exy_corrected)
            disp('================= exy already corrected ! ========================');
            exy_corrected = 1;
        else
            disp('================= exy being corrected here ! =======================');
            exy = -exy;
            exy_corrected = 1;
        end
        % remove bad data points
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
        
        outlier_removed = 1;
        save([strainFileName,'_v73.mat'],'outlier_removed','exy_corrected','-v7.3');
        
        myFile = matfile(strainFileName);
        myFields = who(myFile);
        for ii=1:length(myFields)
            save([strainFileName,'_v73.mat'],myFields{ii},'-append','-v7.3');
        end
    else
        disp('v7.3 file already exist');
    end
    strainFile{iE} = matfile([strainFileName,'_v73.mat']);
end



%% [1] Get data for selected grain and save.
clear data;
ID_target = 639;    % 1144;

for iE = iE_start:iE_stop
    fName_c2t_result = ['WE43_T6_C1_s',num2str(iE),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMap','stru','clusterNumMapCleaned');
    
    % already processed, so: exy_corrected = 1, outlier_removed = 1
    file_name = [dicPath,'\_',num2str(iE),'_v73.mat']; 
    disp(file_name);
    load(file_name,'exx','exy','eyy','sigma');     % Look at exx, but this can be changed in the future.   % ----------------------------------------------------------------------------------
    
    %%%% grain    
    ID_current=ID_target;
    
    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
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
    sigma_local = sigma(indR_min:indR_max, indC_min:indC_max);
    
    clusterNumMapLocal = clusterNumMap(indR_min:indR_max, indC_min:indC_max);
    clusterNumMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
    exx_local(ID_local~=ID_current) = nan;  % The nans are used for making images for the presentation. 0 is OK in processing, because
    exy_local(ID_local~=ID_current) = nan;
    eyy_local(ID_local~=ID_current) = nan;
    sigma_local(ID_local~=ID_current) = nan;
    
    s.iE = iE;
    s.indR_min = indR_min;
    s.indR_max = indR_max;
    s.indC_min = indC_min;
    s.indC_max = indC_max;
    s.boundaryTF_local = boundaryTF_local;
    s.x_local = x_local;
    s.y_local = y_local;
    s.ID_local = ID_local;
    s.ID_current = ID_current;
    s.exx_local = exx_local;
    s.exy_local = exy_local;
    s.eyy_local = eyy_local;
    s.sigma_local = sigma_local;
    s.clusterNumMapLocal = clusterNumMapLocal;
    
    % cleaned
    clusterNumMapCleanedLocal = clusterNumMapCleaned(indR_min:indR_max, indC_min:indC_max);
    clusterNumMapCleanedLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
    s.clusterNumMapCleanedLocal = clusterNumMapCleanedLocal;

    data(iE) = s;

end

grain_data_file = ['WE43_T6_C1_s_all_grain_',num2str(ID_current),'_local_map.mat'];
save(fullfile(working_dir,grain_data_file),'data');




%% Necessary information for cNumMaps etc.
egc = 506391;    % 211441;   % [iE,grainID,iCluster] format
iC_target = mod(egc,10);
egc = (egc-iC_target)/10;
ID_target = mod(egc,10000);
iE_target = (egc-ID_target)/10000;

for iE = iE_start:iE_stop
    fName_c2t_result = ['WE43_T6_C1_s',num2str(iE),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru');
    struCell{iE} = stru;
end

iS = find(arrayfun(@(x) x.gID == ID_target,struCell{2}));    % just find iS from any of the stru
[iE_list, iC_list] = find_tracked_iE_iC_list(struCell, iS, iE_target, iC_target);


%% (1) plot the strain map, and strain distribution, at strain levels [3:5]

for iE = 5

grain_data_file = ['WE43_T6_C1_s_all_grain_',num2str(ID_target),'_local_map.mat'];
load(fullfile(working_dir,grain_data_file),'data');
s = data(iE);

indR_min = s.indR_min;
indR_max = s.indR_max;
indC_min = s.indC_min;
indC_max = s.indC_max;
boundaryTF_local = s.boundaryTF_local;
x_local = s.x_local;
y_local = s.y_local;
ID_current = s.ID_current;
ID_local = s.ID_local;
exx_local = s.exx_local;
exy_local = s.exy_local;
eyy_local = s.eyy_local;
sigma_local = s.sigma_local;
clusterNumMapLocal = s.clusterNumMapLocal;
clusterNumMapCleanedLocal = s.clusterNumMapCleanedLocal;


%% (1.1) plot/adjust strain map here -------------------------------------------------
[f,a,c] = myplot(x_local-x_local(1), y_local-y_local(1), exx_local);  % caxis([-0.1, 0.00]);  % caxis([-0.14, 0.07]);
caxis([-0.1, 0.00]);
title(a,'\epsilon_x_x');
set(a,'fontsize',24,'xTickLabel',{''},'yTickLabel',{''});
axis equal;
imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_exx_scaleA.tif']);
print(fullfile(working_dir,imgName),'-dtiff');   % to parent folder
close(f);



end


%% (2) Illustrate the [process] of how to determine [nCluster] to use.  Here are some comparisons of using different nCluster to try. 
close all;
colors = parula(6);
colorMap = colors([1,5],:);
colorMap = [.98 .98 .98; colorMap];

rng(0);     % adjust this if you don't like the cluster color

maxCluster = 5;

ind = find((ID_local==ID_current)); % ind = find((ID_local==ID_current)&(~isnan(exx_local))&(~isnan(exy_local))&(~isnan(eyy_local)));
exx_t = exx_local(ind);
exy_t = exy_local(ind);
eyy_t = eyy_local(ind);
data_t = [exx_t(:), exy_t(:), eyy_t(:)];

% first, predict centroid
[~,ind_centroid_initial] = max(stru(iS).tSF);
centroid_initial = stru(iS).tStrain(ind_centroid_initial,:);

% ======================= kmeans, determine optimum number of clusters ====================================
% The way it was actually used.
nPoints = 10000;
ind_reduce = ~isnan(sum(data_t,2));
data_reduce = data_t(ind_reduce,:);
reduce_ratio = ceil(size(data_reduce,1)/nPoints);
data_reduce = data_reduce(1:reduce_ratio:end,:);

% To fit here, make number of data points small enough to make the scatter plot :
% data_reduce = data_t(1:step:end,:);

if(~isempty(data_reduce))
    % compare the silhouette, by actually do kmeans on down-sampled samples.
    disp(['ID=',num2str(ID_current)]);
    clear wssd  score_avg  score_cluster_mean  score_cluster_neg_sum  score_cluster_mean_min  score_neg_sum;
    %         score_min = -1*ones(1, maxCluster);
    score_neg_sum = -inf*ones(1, maxCluster);
    nRep = 1;
    c0 = kmeans_pp_init(data_reduce,maxCluster,nRep,centroid_initial);
    
    for nc = 2:2    %maxCluster
        [idx, centroid, sumd] = kmeans(data_reduce, nc, 'Distance','sqeuclidean','MaxIter',1000,'start',c0(1:nc,:,:));   % 'correlation' distance not good.
        sil_score = silhouette(data_reduce,idx);
        
        score_avg(nc) = nanmean(sil_score); % avg score for the condition of nc clusters
                
        % record the vol pct of each cluster, in order to match clusterNum later
        clear cvPctK;
        for ii=1:nc
            cvPctK(ii) = sum(idx(:)==ii)/sum(idx>0);
        end
        
        % the Silhouette figures:

        % (2) with average silhouette
        figure; silhouette(data_reduce,idx);
        set(gca,'fontsize',18,'xlim',[-0.5,1]);
        % also, add the average silhouette value as a vertical line
        hold on; fplot(@(x) (x-score_avg(nc))*2^64,'--r','linewidth',2);
        
        title(['K = ',num2str(nc),', Avg Silhouette = ',num2str(score_avg(nc),'%.3f')],'fontweight','normal');
        
        imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_silhouette_nc_',num2str(nc),'_withAvg.tif']);
        print(fullfile(working_dir,imgName),'-dtiff');   % to parent folder
        
        % (3) scatter plot show the distribution of the clusters
        figure;
        a = axes;
        hold on;
        plot3(a,stru(iS).tStrain(:,1),stru(iS).tStrain(:,2),stru(iS).tStrain(:,3),'.k','markersize',24);
        
         for ii = 1:nc
            xx = data_reduce(idx==ii,1);
            yy = data_reduce(idx==ii,2);
            zz = data_reduce(idx==ii,3);
            %     plot3(xx(1:step:end),yy(1:step:end),zz(1:step:end),'.','color',colors(ii,:));
            scatter3(a,xx,yy,zz,10,'MarkerEdgeColor','none','MarkerFaceColor',colorMap(ii+1,:),'MarkerFaceAlpha',0.33);
        end
        xlabel('\epsilon_x_x');
        ylabel('\epsilon_x_y');
        zlabel('\epsilon_y_y');
        set(a,'fontsize',18,'xlim',[-0.1,0.05],'ylim',[-0.07,0.08],'zlim',[-0.05 0.1]);
        
        view([0 0 1]);

        view(25,20);
        imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_eij_scatter_nc_',num2str(nc),'.tif']);
        print(fullfile(working_dir,imgName),'-dtiff');   % to parent folder
        
        % (4) clusterNumMap if using this nCluster
        [idxx, centroid, sumd] = kmeans(data_t, nc, 'Distance','sqeuclidean','MaxIter',1000,'start',c0(1:nc,:,:)); 
        % reorder idx
        clear cvPct;
        for ii=1:nc
           cvPct(ii) = sum(idxx==ii)/sum(idxx>0);
        end
        [~,c_template] = sort(cvPctK,'descend');
        [~,c_this] = sort(cvPct,'descend');
        for ii=1:nc
           idxx(idxx==(c_this(ii))) = c_template(ii)+100; 
        end        
        idxx = idxx - 100;
        
        cNumMapTemp = zeros(size(ID_local));    % temporary cluster number map  
        cNumMapTemp(ind) = idxx;                 % copy cluster #
        cNumMapTemp(ID_local~=ID_current) = 0;  % remove non-grain data 
        cNums = unique(cNumMapTemp(:));
        cNums(isnan(cNums)) = [];
        n = length(cNums);
        
        % [plot] (temp) cluster ID map
        [f,a,c] = myplot(cNumMapTemp);
        colormap(colorMap);
        caxis(caxis+[-0.5 0.5]);
        c.Ticks = 1:n;
        set(c,'limits',[0.5,n-0.5]);
        title(a,'Cluster ID', 'fontweight','normal');
        set(a,'fontsize',24,'xTickLabel',{''},'yTickLabel',{''});
        axis equal;
        imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_clusterNumMap_with_K_',num2str(nc),'.tif']);
        print(fullfile(working_dir,imgName),'-dtiff');   % to parent folder
        % close(f);
    
    end
    
end

%% [1] cluster ID map, [2] clean up, and cleaned cluster ID map, [3] highlight cluster of interest as red
% plot (temp) cluster ID map
[f,a,c] = myplot(cNumMapTemp);
colormap(colorMap);
caxis([-0.5,2.5]);
set(c,'limits',[0.5,2.5]);
c.Ticks = 1:2;
title(a,' ', 'fontweight','normal');
set(a,'fontsize',24,'xTickLabel',{''},'yTickLabel',{''});
axis equal;
imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_clusterNumMap_with_K_',num2str(nc),'.tif']);
print(fullfile(working_dir,imgName),'-dtiff');   % to parent folder

% cleaned cluster ID map
cNumMapCleaned = one_pass_fill_and_clean(cNumMapTemp, 0.00025);
cNumMapCleaned(ID_local~=ID_current)=0;
[f,a,c] = myplot(x_local, y_local, cNumMapCleaned, (boundaryTF_local));
colormap(colorMap);
caxis([-0.5,2.5]);
set(c,'limits',[0.5,2.5]);
c.Ticks = 1:2;
title(a,' ', 'fontweight','normal');
set(a,'fontsize',24,'xTickLabel',{''},'yTickLabel',{''});
axis equal;
imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_clusterNumMapCleaned_with_K_',num2str(nc),'.tif']);
print(fullfile(working_dir,imgName),'-dtiff');   % to parent folder




%% thinning

clusterNumMapC = cNumMapCleaned;
clusterNumMapC(clusterNumMapC~=iC_target) = 0;
    
clusterNumMapT = double( bwskel(imbinarize(clusterNumMapC),'MinBranchLength',0 * round(min(size(clusterNumMapC))*0.05)) );
t = clusterNumMapT;
t(t==0) = -inf;
myplotm(t,boundaryTF_local)
ca([0, 10]);
title('skeleton of cluster considered');
set(gca,'xticklabel','','yticklabel','','fontsize',18);
title('Skeleton', 'fontweight', 'normal'); title('');
colorbar('off');
print(fullfile(working_dir, 'thinned_cluster_map.tiff'), '-dtiff');

%% hough transform
% (6) Then do hough transform. H = intensity on [rhos, thetas] map
[H,Theta,Rho] = hough(clusterNumMapT,'RhoResolution',1);    % H: d-1 is rho, d-2 is theta.

% (7) Find peaks. Set a [neighborhood size], d_width = 5% map size, d_angle = 5 deg.
maxNumPeaks = 32;
peaks = houghpeaks(H, maxNumPeaks, 'Threshold', 0.3 * max(H(:)), 'NHoodSize',[round_odd(0.05*min(size(clusterNumMapC))),5] );
peakAngles = Theta(peaks(:,2));
peakStrength = H(sub2ind(size(H),peaks(:,1),peaks(:,2)));
disp( table(peakAngles(:),peakStrength(:),'VariableNames',{'PeakAngles','PeakStrength'}) );


% This is just to [illustrate] where the peak is in the hough space
myplot(Theta,Rho,H);
title('Hough Transform of Skeleton','fontweight','normal'); title('');
axis normal;
hold on;
for k = 1:size(peaks,1)
    xy = [Theta(peaks(k,2)), Rho(peaks(k,1))];
    plot3(xy(1),xy(2),max(H(:)),'s','LineWidth',((maxNumPeaks+1-k)/maxNumPeaks)*4,'Color','k');
end
xlabel('\theta, degrees'); ylabel('\rho');
set(gca,'fontsize',24,'xTick',[-90:45:90]);
print(fullfile(working_dir, 'hough_transform.tiff'), '-dtiff');




%% [Show the New method] cluster to twin variant, for iE_select, iS, iC_select  
% The main code was in function: ..._confirmed_twin_to_variant_()
% This is done after twin identification. So, we need to use
% 'clusterNumberMapCleaned' + 'struCell{iE}(iS).cTrueTwin'

fName_c2t_result = ['WE43_T6_C1_s',num2str(iE),'_cluster_to_twin_result.mat'];
load([saveDataPath,fName_c2t_result],'clusterNumMap','clusterNumMapCleaned');
load('D:\p\m\DIC_Analysis\temp_results\WE43_T6_C1_new_variant_map_20200401.mat', 'struCell');

stressTensor = [-1,0,0; 0,0,0; 0,0,0];
ID_current = struCell{iE}(iS).gID;
ind = find(gID==ID_current);

euler = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
[abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [0,0,0], [0,0,0], stressTensor, 'Mg', 'twin');

[ssa, c_a, nss, ntwin, ssGroup] = define_SS('Mg','twin');
traceDir = abs_schmid_factor(nss+1:nss+ntwin,3);    % angle x-to-y

nNeighbors = gNNeighbors(ind);
ID_neighbors = gNeighbors(ind, 1:nNeighbors);

ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);

% Make it one data point wider on each side
indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);

% (Step-1) Crop local maps
ID_local = ID(indR_min:indR_max, indC_min:indC_max);
X_local = X(indR_min:indR_max, indC_min:indC_max);
Y_local = Y(indR_min:indR_max, indC_min:indC_max);
uniqueBoundary_local = uniqueBoundary(indR_min:indR_max, indC_min:indC_max);
uniqueBoundary_local((floor(uniqueBoundary_local/10000)~=ID_current)&(mod(uniqueBoundary_local,10000)~=ID_current)) = 0;    % leave only associated with this grain.
boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);

clusterNumMap_local = clusterNumMapCleaned(indR_min:indR_max, indC_min:indC_max);
clusterNumMap_local(ID_local~=ID_current) = 0;
% Find active system, if any, using cTrueTwin/tGb field
activeTS = sum(struCell{iE}(iS).cTrueTwin,1)>0;

grains_variant_map = zeros(size(ID_local));
% for each cluster, depending on how many variants it has
for iC = 1 % 1:size(struCell{iE}(iS).cTrueTwin,1) % iC_target 
    
    nTwins = sum(struCell{iE}(iS).cTrueTwin(iC,:));
    twinClusters = iC;
    
    if nTwins==0
        % not a twinned cluster
        tnMap = zeros(size(clusterNumMap_local));
    elseif nTwins==1
        % assign this cluster to a single variant
        variant_num = find(struCell{iE}(iS).cTrueTwin(iC,:));
        tnMap = zeros(size(clusterNumMap_local));
        tnMap(clusterNumMap_local==iC) = variant_num;
    elseif nTwins>1
        clear csl;
        
        % use red color to indicate the cluster being analyzed
        map_this_cluster = double(clusterNumMap_local==iC);
        map_this_cluster(map_this_cluster==0) = nan;
        [f,a,c] = myplot(X_local,Y_local,map_this_cluster,boundaryTF_local);
        temp_cmap = [1,1,1; 1,0,0];
        colormap(temp_cmap);
        alpha(0.4);
        colorbar off; grid off;
        set(gca,'xTick',[],'yTick',[]);
        title('');
        print(fullfile(working_dir, 'cluster_analyzed.tiff'),'-dtiff');
        
        % need to divide
        for iTwin = 1:6
%             if activeTS(iTwin) == 1 % [option-2]
            if struCell{iE}(iS).cTrueTwin(iC,iTwin)==1
                % (Step-2) Rotate maps
                ID_r = imrotate(ID_local, traceDir(iTwin), 'nearest', 'loose');
                X_r = imrotate(X_local, traceDir(iTwin), 'nearest', 'loose');
                Y_r = imrotate(Y_local, traceDir(iTwin), 'nearest', 'loose');
                uniqueBoundary_r = imrotate(uniqueBoundary_local, traceDir(iTwin), 'nearest', 'loose');
                uniqueBoundary_r = imdilate(uniqueBoundary_r,ones(3));
                boundaryTF_local_r = imrotate(boundaryTF_local, traceDir(iTwin), 'nearest', 'loose');
                boundaryTF_local_r = imdilate( boundaryTF_local_r,ones(3));
                clusterNumMap_r = imrotate(clusterNumMap_local, traceDir(iTwin), 'nearest', 'loose');
                vMap_r = ismember(clusterNumMap_r, twinClusters); % try this to segment clusterNumMap into variantMap(or trueTwinMap)
                
                % myplot(vMap_r, uniqueBoundary_r);
                
                [nr,nc] = size(vMap_r);
                gbLabelMap = zeros(nr,nc);    % to store assigned gb_label
                gbNumXY_intersect = [];     % [gbNum, Xpos, Ypos]  ----------------------------------> this could be recorded in struCell.
                cslMap = zeros(nr,nc);     % a map recording Connected Segment Length (CSL)   ----------> this might be helpful do determine variant number, if starting point is clusterNumMap.
                gbLR = zeros(nr,2);     % store each rows two possible gbs
                
                % (Step-3)
                for ir = 1:nr
                    if any(vMap_r(ir,:))
                        icL_back = find(uniqueBoundary_r(ir,:),1,'first');
                        gbL = uniqueBoundary_r(ir,icL_back);
                        icR_back = find(uniqueBoundary_r(ir,:),1,'last');
                        gbR = uniqueBoundary_r(ir,icR_back);
                        gbLR(ir,:) = [gbL, gbR];
                        
                        % (instert Step-4) determine if gbL/R can be considered as an intersecting gb. -----------------------------------
                        length_cr = round(min(30, (icR_back - icL_back)/2));
                        num_cr = round(length_cr * 0.7);
                        if sum(vMap_r(ir,icL_back:icL_back+length_cr))>num_cr
                            gbNumXY_intersect = [gbNumXY_intersect; gbL, X_r(ir,icL_back), Y_r(ir,icL_back)];
                        end
                        if sum(vMap_r(ir,icR_back-length_cr:icR_back))>num_cr
                            gbNumXY_intersect = [gbNumXY_intersect; gbR, X_r(ir,icR_back), Y_r(ir,icR_back)];
                        end
                        % end of (Step-4). Alternatively, read from manual label for data analysis. ---------------------------------------
                        
                        icL_front = find(vMap_r(ir,:),1,'first');
                        icR_front = find(vMap_r(ir,:),1,'last');
                        % will be false, if either is empty
                        while (icL_front<=icR_front)
                            csl_length = 0;
                            if (icL_front-icL_back)<=(icR_back-icR_front)
                                % search for connected segments from left to right
                                while(vMap_r(ir,icL_front))
                                    gbLabelMap(ir,icL_front) = 1;  % prepare to assign label
                                    vMap_r(ir,icL_front) = 0; % make element on variant map 0
                                    icL_front = icL_front + 1;   % move pointer forward to the right
                                    csl_length = csl_length + 1;
                                end
                                cslMap(gbLabelMap==1) = csl_length;
                                gbLabelMap(gbLabelMap==1) = gbL;
                                icL_back = icL_front - 1;    % assign left side back
                                icL_front = find(vMap_r(ir,:),1,'first'); % search left side front again
                            else
                                while (vMap_r(ir,icR_front))
                                    gbLabelMap(ir,icR_front) = 1;
                                    vMap_r(ir,icR_front) = 0;
                                    icR_front = icR_front - 1;
                                    csl_length = csl_length + 1;
                                end
                                cslMap(gbLabelMap==1) = csl_length;
                                gbLabelMap(gbLabelMap==1) = gbR;
                                icR_back = icR_front + 1;
                                icR_front = find(vMap_r(ir,:),1,'last');
                            end
                        end
                    end % end of if any(variant_r(ir,:))
                end % end of for ir=1:nr
                
                % (Step-5) go back to clean.  Sometimes, no intersection is determined ...
                cleanTF = 1;
                if (cleanTF)&&(~isempty(gbNumXY_intersect))
                    gbList = unique(gbNumXY_intersect(:,1));
                    for ir=1:nr
                        tf = ismember(gbLR(ir,:),gbList);
                        if sum(tf)==1
                            gbOK = gbLR(ir, tf);
                            gbKO = gbLR(ir, ~tf);
                            ind = gbLabelMap(ir,:) == gbKO;
                            gbLabelMap(ir,ind) = gbOK;
                        elseif sum(ismember(gbLR(ir,:),gbList))==0
                            gbLabelMap(ir,:) = 0;
                        end
                    end
                end
                
                cslMap_t = cslMap;
                cslMap_t(cslMap_t==0) = nan;
                myplotm(cslMap_t,'TF',boundaryTF_local_r)
                ca([0 370]);
                set(gca,'fontsize',24);xlabel('x');ylabel('y');
                title('');

                
                temp = cslMap;
                % rotate back, need to crop again.
                temp = imrotate(cslMap,-traceDir(iTwin), 'nearest', 'loose');
                ID_back = imrotate(ID_r,-traceDir(iTwin), 'nearest', 'loose');
                ind_back = ismember(ID_back, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
                
                % Need to crop a region of the original size.
                img1_template = (ID_local==ID_current);
                img2_signal = (ID_back==ID_current);
                [yOffSet, xOffSet] = normxcorr2A_register(img1_template, img2_signal, [0 0 0 0], [0 0 0 0], 0);
                indC_back_min = 1 + xOffSet;
                indC_back_max = indC_back_min + indC_max - indC_min;
                indR_back_min = 1 + yOffSet;
                indR_back_max = indR_back_min + indR_max - indR_min;
                
                csl(:,:,iTwin) = temp(indR_back_min:indR_back_max, indC_back_min:indC_back_max);
                
                temp = csl(:,:,iTwin);
                temp(clusterNumMap_local~=iC) = nan;
                myplotm(temp, imdilate(boundaryTF_local,ones(3)));
                caxis([0 370]);
                set(gca,'xTick',[],'yTick',[],'fontsize',18);
                
            end % end of if(activeTS(iTwin)==1)
            
        end % end of for iTwin=1:6
        
        [nr,nc,np] = size(csl);
        tnMap = zeros(nr,nc);   % variant_num_map of this cluster
        for ir=1:nr
            for ic=1:nc
                [maxV,temp] = max(csl(ir,ic,:));
                if maxV>0
                    tnMap(ir,ic) = temp;
                end
            end
        end
        
    end % end of elseif sum(struCell{iE}(iS).cTrueTwin(iC,:))>1
    
    % Here need to prevent elements belonging to other clusters having a tnMap value (caused by rotation).
    tnMap(~ismember(clusterNumMap_local,twinClusters)) = 0;
    
    tnMap_t = tnMap;
    tnMap_t(tnMap_t==0) = -inf;
    [f,a,c] = myplot(tnMap_t,boundaryTF_local);
    colormap(colors([1,5],:));
    
    set(c,'limits',[1-1.5, 4+1.5],'Ticks',[1, 4],'TickLabels',{'1','4'});
    set(gca,'xticklabel','','yticklabel','');
    title('Variant ID, uncleaned','fontweight','normal');title('');
    set(gca,'xTick',[],'yTick',[],'fontsize',24);
    print(fullfile(working_dir,'raw_variantMap.tiff'),'-dtiff');

    
    clusters_variant_map{iC} = tnMap;
    grains_variant_map = grains_variant_map + clusters_variant_map{iC};
end % end of for iC = 1:num_of_clusters

% do clean up here, after combining all the clusters.  
% Note: if not all clusters are selected, the result can be different from that actually analyzed for the data.  
grains_variant_map  = one_pass_fill(grains_variant_map);

grains_variant_map(grains_variant_map==0) = -inf;   % -> make background white
[f,a,c] = myplot(grains_variant_map,boundaryTF_local);
colormap(colors([1,5],:));

set(c,'limits',[1-1.5, 4+1.5],'Ticks',[1, 4],'TickLabels',{'1','4'});
set(gca,'xticklabel','','yticklabel','');
title('Variant ID, cleaned','fontweight','normal');title('');
set(gca,'xTick',[],'yTick',[],'fontsize',24);
print(fullfile(working_dir,'cleaned_variantMap.tiff'),'-dtiff');

close all;




