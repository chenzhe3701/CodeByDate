
%% set up
clear; close all; clc;
addChenFunction;
data_dir = 'E:\Mg4Al_S1_insitu\Analysis_by_Matlab';
variant_data = 'E:\Mg4Al_S1_insitu\Analysis_by_Matlab\previous result recovered\20200325_0506_new_variant_map_Mg4Al_S1.mat';
dic_dir = 'E:\Mg4Al_S1_insitu\SEM Data\stitched_DIC';
sample_name = 'Mg4Al_S1';

output_dir = 'E:\zhec umich Drive\0_temp_output\Mg4Al_S1 analysis for paper';
mkdir(output_dir);

strains = [0, -0.0012, -0.0117, -0.0186, -0.0172, -0.0124, 0.0006]; % iE=0:6

%% load data
d = matfile(fullfile(data_dir, 'Mg4Al_S1_EbsdToSemForTraceAnalysis.mat'));
ID = d.ID;
X = d.X;
Y = d.Y;
boundaryTF = d.boundaryTF;
gID = d.gID;
gPhi1 = d.gPhi1;
gPhi = d.gPhi;
gPhi2 = d.gPhi2;
gDiameter = d.gDiameter;

d2 = matfile(variant_data);
variant_cell = d2.variantMapCleanedCell;
struCell = d2.struCell;

clear cluster_number_cell;
for iE = 1:6
    d3 = matfile(fullfile(data_dir, ['Mg4Al_S1_s',num2str(iE),'_cluster_to_twin_result.mat']));
    cluster_number_cell{iE} = d3.clusterNumMapCleaned;
end

%% convert iE to iB for consistency
for iE=6:-1:0
    iB = iE + 1;
    if iE==0
        variant_cell{iB} = zeros(size(X));
        cluster_number_cell{iB} = zeros(size(X));
    else
        variant_cell{iB} = variant_cell{iE};
        cluster_number_cell{iB} = cluster_number_cell{iE};
    end
end

%% generate twin map: segment variant map, and assign unique twin grain ID 
for iE = 0:6
    iB = iE + 1;
    twin_map = variant_cell{iB};
    cluster_map = cluster_number_cell{iB};
    
    temp_map = ID*100 + cluster_map*10 + twin_map;   % unique label for each connected area
    
    twin_grain_ID = one_pass_label(temp_map);   % one-pass-label twin grains
    twin_grain_ID(twin_map<1) = 0;       % make untwinned area 0
    
    % check how many grains
%     uniqueIDs = nan_unique(twin_grain_ID(:));
%     uniqueIDs(isnan(uniqueIDs)) = [];
%     length(uniqueIDs)  
        
    twin_grain_cell{iB} = twin_grain_ID;
end

%% ==> consider reducing size of map by 4
reduce_ratio = 4;
X = X(1:reduce_ratio:end, 1:reduce_ratio:end);
Y = Y(1:reduce_ratio:end, 1:reduce_ratio:end);
ID = ID(1:reduce_ratio:end, 1:reduce_ratio:end);
boundaryTF = find_one_boundary_from_ID_matrix(ID);

for iE = 0:6
   iB = iE + 1;
   variant_cell{iB} = variant_cell{iB}(1:reduce_ratio:end, 1:reduce_ratio:end);
   cluster_number_cell{iB} = cluster_number_cell{iB}(1:reduce_ratio:end, 1:reduce_ratio:end);
   twin_grain_cell{iB} = twin_grain_cell{iB}(1:reduce_ratio:end, 1:reduce_ratio:end);
end

% plot
for iE = 0:6
   iB = iE + 1;      
   myplot(X,Y,rem(twin_grain_cell{iB},5) + double(twin_grain_cell{iB}>0), boundaryTF);
   set(gca,'xTickLabel',[],'yTickLabel',[]);
   print(fullfile(output_dir,['twin_grain_map_iE=',num2str(iE),'.tiff']),'-dtiff');
   close;
end

save(fullfile(output_dir, 'temp_ws_1.mat'));
%% Categorize twin data, on pixel level, and on grain level
% pixel level, at each iE, there are:
% (p1) current twin (directly from twin map)
% (p2) past-or-present twin (obtained by accummulating twin map up to this iE) 
% based on these two basic maps, we can get
% (p3, yellow) fresh-twin = currently twinned - ever twinned before = p1 - p2(iE-1);  
% (p4, blue) recurring-twin = currently twinned & ever twinned bofore = p1 & p2(iE-1); 
% (p5, gray) de-twin = ever twinned before - currently twinned = p2(iE-1) - p1;
% p1 = p3 + p4
% p2 = p3 + p4 + p5
% 
% So, we have these pixel level maps:
% [1] p1 map, current twin map;
% [2] p2 map, ever twinned map; (p1 is a subset of p2)
% [3] p3 + p4 + p5 map, fre/recurring/de-twin map; 
p0 = 0;
p1 = 1;
p2 = 2;
p3 = 3;
p4 = 4;
p5 = 5;

% For SEM-DIC sample, maps of different iEs are already overlapped, so no need to define again a valid region   
for iE = 0:6
    iB = iE + 1;
    current_twinTF_cell{iB} = double(variant_cell{iB}>0);
end

for iE = 0:6
    iB = iE + 1;
    if iE==0
        past_present_twin_cell{iB} = zeros(size(ID));    % max-twin area up to iE
        fresh_twin_cell{iB} = zeros(size(ID));    % new-twin at iE. max_twin(iE-1) + new_twin(iE) = max_twin(iE)
        recurring_twin_cell{iB} = zeros(size(ID));  % re_twin(iE) + new_twinned(iE) = twin(iE)
        de_twin_cell{iB} = zeros(size(ID));  % max_twin(iE) - de_twin(iE) = twin(iE)
    else
        past_present_twin_cell{iB} = double(past_present_twin_cell{iB-1} | current_twinTF_cell{iB}>0);
        fresh_twin_cell{iB} = double(current_twinTF_cell{iB}>0 & past_present_twin_cell{iB-1}==0);
        recurring_twin_cell{iB} = double(current_twinTF_cell{iB}>0 & past_present_twin_cell{iB-1}>0);
        de_twin_cell{iB} = double(current_twinTF_cell{iB}==0 & past_present_twin_cell{iB-1}>0);
    end
end
% new/re/de-twin map: p3=new-twin, p4=re-twin, p5=de-twin
for iE = 0:6    
    iB = iE + 1;
    
    frd_twin_cell{iB} = fresh_twin_cell{iB}*p3 + recurring_twin_cell{iB}*p4 + de_twin_cell{iB}*p5;
    
    % Plot to illustrate
    [f,a,c] = myplot(X,Y,frd_twin_cell{iB}, boundaryTF);
        
    colors = parula(16);
    cmap = [1 1 1;      % backgrond, value < 2
        colors(14,:);   % yellowis, new twin, value = 3
        colors(2,:);    % blueish, retwin, value=4
        .5, .5, .5];    % gray, detwin, value=5
    caxis([1.5, 5.5]);
    colormap(cmap);
    
    set(gca,'xTickLabel',[],'yTickLabel',[]);
    set(c,'limits',[2.5,5.5], 'Ticks',[3,4,5], 'TickLabels',{'fresh-twin', 'recurring-twin', 'de-twin'}); 
    title(['pixel summary, load step = ',num2str(iE), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','norma');
    print(fullfile(output_dir,['frd twinned iE=',num2str(iE),'.tiff']),'-dtiff');
    
    close all;
end


%% Generating grain level map
% current twin grain:
g1 = 1; % (g1) fresh twin
g2 = 2; % (g2) de-twin then re-twin
g3 = 3; % (g3) evolving current twin
% ever twin grain
g4 = 4; % (g4) evolving past-or-present twin
g5 = 5; % (g5) completely de-twin

[nR,nC] = size(ID);

for iE = 0:6
    iB = iE + 1;
    
    current_twin_grain_label = zeros(nR,nC);    
    if iE==0
        past_or_present_twin_grain_ID = zeros(nR,nC);
        past_or_present_twin_grain_label = zeros(nR,nC);
    else
        past_or_present_twin_grain_ID = past_or_present_twin_grain_ID_cell{iB-1};
        past_or_present_twin_grain_label = past_or_present_twin_grain_label_cell{iB-1};
    end
    ID_assign = max(past_or_present_twin_grain_ID(:)) + 1;  % initial ID to be assigned to ever twin grain 
    
    ID_c = twin_grain_cell{iB};  
    
    % [step-1] loop check every 'current twin grain'
    twin_grain_list = nan_unique(ID_c(:));
    twin_grain_list(isnan(twin_grain_list)) = [];   % remove nan
    twin_grain_list(twin_grain_list==0) = [];       % remove 0 (if any)
    
    disp(['iE=',num2str(iE),', # of grains=',num2str(length(twin_grain_list))]);
    pause(10);
    for ii=1:length(twin_grain_list)
       ID_c_current = twin_grain_list(ii);
             
       % index of ID_c_current
       inds = ismember(ID_c, ID_c_current);
       
       % If this grain is a twin grain (contains twin pixel)
       grain_twinned_TF = any(current_twinTF_cell{iB}(inds));
       if grain_twinned_TF
           % determine twin map labels: g1 = new twin grain, g2 = mixed
           labels = frd_twin_cell{iB-1}(inds);
           if sum(labels==p0)/numel(labels) >0.9    % intersect(labels, [p0,p3,p4,p5]) == p0
               % the grain is completely new twin (g1/new twin), as it does not overlap with any previously twinned pixel  
               current_twin_grain_label(inds) = g1;
           elseif intersect(labels, [p3,p4,p5]) == p5
               % if most are p5, we can allow a little bit p0 pixels
               current_twin_grain_label(inds) = g2;
           else
               % evolving grain
               current_twin_grain_label(inds) = g3;
           end
           
           % current twin grain will contribute to the 'ever twin grain map' 
           ids = unique(past_or_present_twin_grain_ID(inds));
           ids(ids==0) = [];
           if isempty(ids)
               % if current twin grain does not overlap with any ever-twinned grain, ADD this grain to the ever twin grain map 
               past_or_present_twin_grain_ID(inds) = ID_assign;
           else
               % if this grain overlap with ever twinned grain, MERGE them. (modify inds to include all grains)  
               inds = ismember(past_or_present_twin_grain_ID, ids) | inds;
               past_or_present_twin_grain_ID(inds) = ID_assign;
           end
           ID_assign = ID_assign + 1;
           
           % this 'merged ever twin grain' contains 'current twin', so at least p4/re-twin, maybe p3/new-twin pixels
           past_or_present_twin_grain_label(inds) = g4;    % 'g4/evolving'
       end
    end
    current_twin_grain_label_cell{iB} = current_twin_grain_label;
    
    % [step-2] loop check every 'ever twin grain', find the completely detwinned grain ===> maybe need to move this to [step-1] 
    pp_twin_grain_list = nan_unique(past_or_present_twin_grain_ID(:));
    pp_twin_grain_list(isnan(pp_twin_grain_list)) = []; % remove nan
    pp_twin_grain_list(pp_twin_grain_list==0) = [];     % remove 0 (if any)  
    
    for ii = 1:length(pp_twin_grain_list)
        ID_current = pp_twin_grain_list(ii);
        inds = ismember(past_or_present_twin_grain_ID, ID_current);        
           
        % make ever_twin_grain_ID and ever_twin_grain_label
        labels = unique(frd_twin_cell{iB}(inds));
        if intersect(labels, [p3,p4,p5]) == p5
            % only contains conpletely detwinned pixels (p5), this ever twin grain should be labeled by (g4)   
            past_or_present_twin_grain_label(inds) = g5;    % 'g5/completely detwin'
        else
            % this currently twinned contains p3/new-twin or p4/re-twin pixels.   
            past_or_present_twin_grain_label(inds) = g4;    % 'g4/mixed'
        end
    end
    
    if iE>0
        % remake IDs
        past_or_present_twin_grain_ID = hungarian_assign_ID_map(past_or_present_twin_grain_ID, past_or_present_twin_grain_ID_cell{iB-1});
    end
    past_or_present_twin_grain_ID_cell{iB} = past_or_present_twin_grain_ID;
    past_or_present_twin_grain_label_cell{iB} = past_or_present_twin_grain_label;
    
end

%% plot summary maps, overlay pixel info + grain info
close all;
for iE = 0:6
    iB = iE + 1;
    
    map = zeros(nR,nC);
    
    inds = past_or_present_twin_grain_label_cell{iB}>0;
    map(inds) = past_or_present_twin_grain_label_cell{iB}(inds);
    
    inds = current_twin_grain_label_cell{iB}>0;
    map(inds) = current_twin_grain_label_cell{iB}(inds);
    
    [f,a,c] = myplot(X,Y, map, boundaryTF);
    
    colors = plasma(25);
    cmap = zeros(6,3);
    cmap(1,:) = [1,1,1];    % 0 = untwinned, background
    cmap(2,:) = colors(24,:);    % g1 = completely new twin, yellowish
    cmap(3,:) = colors(10,:);     % g2 = de-twin then re-twin, purple
    cmap(4,:) = colors(16,:);     % g3 = evoling current twin, pink
    cmap(5,:) = [0.7, 0.7, 0.7];     % g4 = evolving past-or-present twin, light gray
    cmap(6,:) = [0.2, 0.2, 0.2];    % g5 = completely de-twin, dark gray
    
    colormap(cmap);
    caxis([-0.5, 5.5]);
    set(c,'limits',[0.5, 5.5], 'Ticks', 1:5, ...
        'TickLabels', {'g1:fresh','g2:de-twin then re-twin','g3:evolving current','g4:evolving past or present','g5:completely detwin'})
    title(['grain summary, load step = ',num2str(iE), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','norma');
    set(gcf,'position',[80,172,1000,600]);
    set(gca,'xTickLabel',[],'yTickLabel',[]);
    print(fullfile(output_dir,['grain label iE=',num2str(iE),'.tiff']),'-dtiff');
    close;
end

save(fullfile(output_dir, 'temp_ws_2.mat'));

%% speicial make for iE=6, with only 3 types of grain labels (as there is no retwinnig for this sample)
for iE = 4:6
    iB = iE + 1;
    
    map = zeros(nR,nC);
    
    inds = past_or_present_twin_grain_label_cell{iB}>0;
    map(inds) = past_or_present_twin_grain_label_cell{iB}(inds);
    
    inds = current_twin_grain_label_cell{iB}>0;
    map(inds) = current_twin_grain_label_cell{iB}(inds);
    
    [f,a,c] = myplot(X,Y, map, boundaryTF);
    
    colors = plasma(25);
    cmap = zeros(6,3);
    cmap(1,:) = [1,1,1];    % 0 = untwinned, background
    cmap(2,:) = colors(24,:);    % g1 = completely new twin, yellowish
    cmap(3,:) = colors(10,:);     % g2 = de-twin then re-twin, purple
    cmap(4,:) = colors(16,:);     % g3 = evoling current twin, pink
    cmap(5,:) = [0.7, 0.7, 0.7];     % g4 = evolving past-or-present twin, light gray
    cmap(6,:) = [0.2, 0.2, 0.2];    % g5 = completely de-twin, dark gray
    
    colormap(cmap);
    caxis([-0.5, 5.5]);
    set(c,'limits',[2.5, 5.5], 'Ticks', 1:5, ...
        'TickLabels', {'Fresh twin','Re-activated twin','Residual twin','Partially detwin','Completely detwin'});
    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','normal');
    set(gcf,'position',[80,172,1000,600]);
    set(gca,'xTickLabel',[],'yTickLabel',[], 'fontsize',16);
    print(fullfile(output_dir,['detwin summary map iE=',num2str(iE),'.tiff']),'-dtiff');
    close;
end



%% iE=1,2,3 compression, 4,5,6 back to tension
%% Twin statistics summary. What we need at max compression:
% Fig 2(a). Variants twinned/not twinned & pct vs. variant_SF.
% Fig 2(b). Grains twinned/not twinned & pct vs. max_basal_SF.
% Fig 2(c)/(d). Max_basal_SF vs max_twin_SF for twinned/not twinned
% Fig 3. Boxplot twin variant area fraction vs variant_SF
% ==> try same for detwin
%
% Fig 10. M' for active variants, all possible variants, distribution ratio
% 
% At each strain level
% Fig 9. Mean effective strain vs. distance_to_gb, for each gb_zone 

% Summary. Based on script for paper part 2. But should also consider
% future use for in-situ EBSD data.

variableNames = {'iE','ID','max_basal_SF','max_twin_SF','activeTwin_SF','gDia','twinnedTF','tAF'};
T = cell2table(cell(0,length(variableNames)));
T.Properties.VariableNames = variableNames;

variableNames2 = {'iE','ID','iTwin','variant','variant_SF','vActiveTF','vPct','max_basal_SF','max_twin_SF'};
T2 = cell2table(cell(0,length(variableNames2)));
T2.Properties.VariableNames = variableNames2;

stressTensor = [-1 0 0; 0 0 0; 0 0 0];
% um_per_dp = 60/4096*5;

gList = nan_unique(ID(:));  % gIDs in area of interest
gList(gList==0) = [];
gList(isnan(gList)) = [];

for iE = 1:6
    iB = iE + 1;
    twin_map = variant_cell{iB};
    for ii = 1:length(gList)
        
        ID_current = gList(ii);
        
        inds = ID==ID_current;
        
        % Find area fraction
        vol_grain = sum(inds(:));
        variant_ids = twin_map(inds);
        vol_variants = arrayfun(@(x) sum(ismember(variant_ids, x)), 1:6);
        tAF = sum(vol_variants)/sum(vol_grain);
        
        % If variant is twinned at this strain level.
        activeTS = vol_variants>0;
        
        ind_g = find(ID_current==gID);
        gd = gDiameter(ind_g);
        
        euler_current = [gPhi1(ind_g),gPhi(ind_g),gPhi2(ind_g)];
        [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler_current, [0 0 0], [0 0 0], stressTensor, 'Mg', 'twin');
        max_basal_SF = max(abs_schmid_factor(1:3,2));
        twin_SFs = abs_schmid_factor(19:24,2);
        max_twin_SF = max(twin_SFs);
        
        [~,SF_order]=ismember(twin_SFs,sort(twin_SFs,'descend'));   % find index of ranking of element in array
        
        if any(activeTS)
            activeTwin_SF = max(twin_SFs(activeTS));
            twinnedTF = true;
        else
            activeTwin_SF = nan;
            twinnedTF = false;
        end
        
        T = [T; {iE, ID_current, max_basal_SF, max_twin_SF, activeTwin_SF, gd, twinnedTF, tAF}];
        
        for iTwin=1:length(activeTS)
            if activeTS(iTwin)>0
                vActiveTF = true;
            else
                vActiveTF = false;
            end
            vi = SF_order(iTwin);
            vSF = twin_SFs(iTwin);
            vPct = vol_variants(iTwin)/vol_grain;
            T2 = [T2; {iE, ID_current, iTwin, vi, vSF, vActiveTF, vPct, max_basal_SF, max_twin_SF}];
        end
        
        if rem(ii,10) == 1
           str = sprintf('iE=%d, ID=%d', iE, ID_current);
           disp(str);
        end
    end
    
end

%% [plot] histogram showing distribution of variant SF
ind = T2.iE==1;
edges = -0.5:0.05:0.5;
figure;
histogram(T2.variant_SF(ind), edges);
xlabel('Variant Schmid Factor');
ylabel('Counts');
set(gca, 'xTick',-0.5:0.1:0.5, 'fontsize',16);
print(fullfile(output_dir, 'twin variant SF.tiff'), '-dtiff');

%% [paper Fig 9]. Counts of variants twinned & not-twinned vs. variant_SF
edges = -0.5:0.05:0.5;
xstr = [];
for ii=1:length(edges)-1
    xstr{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
end
for iE=1:6
    iB = iE + 1;
    ind = (T2.iE==iE)&(T2.vActiveTF == 1);
    ind2 = (T2.iE==iE)&(T2.vActiveTF == 0);
    [N_t,~] = histcounts(T2.variant_SF(ind), edges);
    [N_nt,~] = histcounts(T2.variant_SF(ind2), edges);
    d_int = (edges(3)-edges(2))/2;
    xpos = [edges(2)-d_int, edges(2:end-1)+d_int];
    
    pct(iE,:) = N_t./(N_t+N_nt);
    
    figure; hold on; disableDefaultInteractivity(gca);
    hbar = bar(xpos, [N_nt(:), N_t(:)], 1, 'stacked');
    % set(gca,'fontsize',12, 'XTick',-0.5:0.1:0.5,'ylim',[0 300]);
    set(gca, 'fontsize',12, 'ylim',[0 300], 'XTick',edges(1:end-1)+d_int, 'xTickLabels',xstr, 'xTickLabelRotation',45, 'xlim',[0 0.5]);
    xlabel('Twin Variant Schmid Factor');
    ylabel('Counts');
    
    yyaxis right;
    set(gca, 'ycolor', 'k','fontsize',16, 'ylim',[0 80]);
    ylabel('Percent (%)');
    plot(-0.475:0.05:0.475, N_t./(N_t+N_nt) * 100,'-ko','linewidth',1.5);

    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','norma');
    title(' ');
    legend({'Variants not twinned', 'Variants twinned','Percent of variants twinned'},'Location','northwest');
    
    hbar(1).FaceColor = [0 0 1];
    hbar(2).FaceColor = [1 0 0];
           
    print(fullfile(output_dir,['fig 9 twin nontwin variant vs SF iE=',num2str(iE),'.tiff']),'-dtiff');

    disp('row1: # twinned, row2: # not twinned, row3: pct');
    clear tt;
    tt = [N_t; N_nt; N_t./(N_t+N_nt)];
    disp(array2table(tt));
    
    close;
end

%% [plot combine of Fig 9] pct variants twinned for all iEs
close all;
colors = inferno(6);
markers = {'-o','-d','-s','-.o','-.d','-.s'};

xstr = [];
for ii=1:length(edges)-1
    xstr{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
end

figure; hold on;
for iE = 1:6
    plot(edges(1:end-1)+0.025, 100*pct(iE,:), markers{iE}, 'color', colors(iE,:), 'linewidth', 1.5);
end
xlabel('Twin Variant Schmid Factor');
ylabel('Percent Twinned (%)');
set(gca, 'fontsize', 16, 'xTick',edges(1:2:end), 'ylim', [0,80]);
legend({'Load step 1','Load step 2','Load step 3','Load step 4','Load step 5','Load step 6'},'location','northwest');

print(fullfile(output_dir,'pct variants twinned vs SF all iEs.tiff'), '-dtiff');

%% [paper Fig 10]. Counts grains twinned & not twinned vs. max_basal_SF
edges = 0:0.05:0.5;
for iE=1:6
    iB = iE + 1;
    ind1 = (T.iE==iE)&(T.twinnedTF==0);
    N_nt = histcounts(T.max_basal_SF(ind1), edges);
    ind2 = (T.iE==iE)&(T.twinnedTF==1);
    N_t = histcounts(T.max_basal_SF(ind2), edges);
    d_int = (edges(3)-edges(2))/2;
    xpos = [edges(2)-d_int, edges(2:end-1)+d_int];
    xstr = [];
    for ii=1:length(edges)-1
        xstr{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
    end
    figure;disableDefaultInteractivity(gca);
    hbar = bar(xpos, [N_nt(:),N_t(:)], 1, 'stacked');
    set(gca,'ylim',[0 100], 'xlim', [0 0.5]);
    xlabel('Maximum Basal Schmid Factor');
    ylabel('Counts');
    
    yyaxis right;
    set(gca,'ycolor','k','XTick',edges(1:end-1)+d_int,'xTickLabels',xstr,'xTickLabelRotation',45);
    plot(xpos, N_t./(N_t+N_nt) * 100,'-ko','linewidth',1.5);
    ylabel('Percent (%)');
    legend('Grains not twinned','Grains twinned','Percent of grains twinned','location','northwest');
    set(gca,'fontsize',12,'ylim',[0 149],'fontsize',16);
    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','norma');
    title('');
    hbar(1).FaceColor = [0 0 1];
    hbar(2).FaceColor = [1 0 0];
        
    print(fullfile(output_dir,['fig 10 twin nontwin grain vs basal SF iE=',num2str(iE),'.tiff']),'-dtiff');
    
    disp(['iE = ',num2str(iE)]);
    tt = [N_t; N_nt; N_t./(N_t+N_nt)];
    disp(array2table(tt'));

end
close all;
%% [Fig 11]
for iE = 1:6
    iB = iE + 1;
    figure;disableDefaultInteractivity(gca); hold on;
    ind = (T.iE==iE)&(T.twinnedTF==0);
    plot(T.max_twin_SF(ind), T.max_basal_SF(ind),'.','markersize',12,'color','b'); % '#0072BD'
    xlabel('Max Twin Schmid Factor');
    ylabel('Max Basal Schmid Factor');
    set(gca,'fontsize',16,'xTick',-0.5:0.1:0.5,'xlim',[-0.5,0.5],'ylim',[0 0.5]);
    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','normal');
    legend('Not twinned','location','southwest');
    
    % title('')

    figure;disableDefaultInteractivity(gca); hold on;
    ind = (T.iE==iE)&(T.twinnedTF==1);
    plot(T.max_twin_SF(ind), T.max_basal_SF(ind),'.','markersize',12,'color','r');    % '#D95319'
    xlabel('Max Twin Schmid Factor');
    ylabel('Max Basal Schmid Factor');
    set(gca,'fontsize',16,'xTick',-0.5:0.25:0.5,'xlim',[-0.5,0.5],'ylim',[0 0.5]);
    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','normal');
    legend('Twinned','location','southwest');
    
    % title('');
    
    figure;disableDefaultInteractivity(gca); hold on;
    ind = (T.iE==iE)&(T.twinnedTF==1);
    plot(T.max_twin_SF(ind), T.max_basal_SF(ind),'.','markersize',12,'color','r');    % '#D95319'
    ind = (T.iE==iE)&(T.twinnedTF==0);
    plot(T.max_twin_SF(ind), T.max_basal_SF(ind),'.','markersize',12,'color','b'); % '#0072BD'
    xlabel('Max Twin Schmid Factor');
    ylabel('Max Basal Schmid Factor');
    set(gca,'fontsize',16,'xTick',-0.5:0.1:0.5,'xlim',[-0.5,0.5],'ylim',[0 0.5]);
    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','normal');
    title(' ');
    legend('Twinned', 'Not twinned', 'location','southwest');
    
    set(gca,'xlim', [ 0, 0.5]); % limit x positive
    axis square;
    
    % title('');
    
    print(fullfile(output_dir,['fig 11 twin nontwin grain vs max SF iE=',num2str(iE),'.tiff']),'-dtiff');
    
    close all;
end

%% [Fig 12]. variant_area_fraction_in_grain vs. variant_SF
for iE = 1:6
    iB = iE + 1;
    ind = (T2.iE==iE)&(T2.vActiveTF==1);
    t2 = T2(ind,:);
    edges = 0:0.05:0.5;
    vx = t2.variant_SF;
    vy = t2.vPct;  
    
    gv = discretize(vx, edges);
    nGroups = length(edges)-1;
    clear labels;
    for ii = 1:length(edges)-1
        labels{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
    end
    figure;disableDefaultInteractivity(gca);
    % boxplot([vy; nan*ones(nGroups, 1)], [gv; (1:nGroups)'],'notch','on');   % prevent having empty group
    boxplot([vy; nan*ones(nGroups, 1)], [gv; (1:nGroups)']);   % prevent having empty group
    
    xlabel('Twin Variant Schmid Factor');
    ylabel('Twin Variant Area Fraction in Grain');
    set(gca,'xticklabels',labels,'xticklabelrotation',45,'fontsize',15);
    set(gca,'ylim',[-0.02 0.62]);
    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','normal');
    title([' ']);
    axis square;
    print(fullfile(output_dir,['fig 12 variant fraction vs SF iE=',num2str(iE),'.tiff']),'-dtiff');
end
close all;

%% summarize/illustrate grain size
close all;
uniqueIDs = unique(ID(:));
for iE = [] %2:6
    iB = iE + 1;
    map_t = zeros(size(ID));                
    
    count_1 = 0;
    size_1 = [];
    count_2 = 0;
    size_2 = [];
    
    for ii=1:length(uniqueIDs)
        ID_current = uniqueIDs(ii);
        inds = ismember(ID, ID_current);
        children = unique(twin_grain_cell{iB}(inds));
        children(children==0) = [];
        children(isnan(children)) = [];
        
        if ~isempty(children)
            parent_size = sum(inds(:));
            for jj = 1:length(children)
                ID_c = children(jj);
                ind_c = ismember(twin_grain_cell{iB}, ID_c);
                child_size = sum(ind_c(:));
                
                % map_t(ind_c) = child_size/parent_size;

                if child_size > 20 %/parent_size > 0.002
                    map_t(ind_c) = 1;
                    count_1 = count_1 + 1;
                    size_1 = [size_1, child_size];
                else
                    map_t(ind_c) = 2;
                    count_2 = count_2 + 1;
                    size_2 = [size_2, child_size];
                end
            end
        end
    end
    myplot(map_t, boundaryTF);
    caxis([0 2]);
    count_1
    count_2
%     size_1 = sort(size_1);
%     size_2 = sort(size_2);
end

%% (A2) variant detwin area (as pct of twin area OR as pct of grain area?) vs. SF 
% This can use pixel map, just compare iE=3 (max) vs iE=6 (min)

variableNames3 = {'iE_ref','iE','ID','iTwin','variant_SF','vActiveTF','A_grain','A_ref','A_def','A_diff'};
T3 = cell2table(cell(0,length(variableNames3)));
T3.Properties.VariableNames = variableNames3;

for iE_ref = 3
    iB_ref = iE_ref + 1;
    twin_map_ref = variant_cell{iB_ref};
    
    for iE = 4:6
        iB = iE + 1;
        twin_map_def = variant_cell{iB};
        
        for ii = 1:length(gList)
            ID_current = gList(ii);
            
            ind_g = find(ID_current==gID);
            gd = gDiameter(ind_g);
            
            euler_current = [gPhi1(ind_g),gPhi(ind_g),gPhi2(ind_g)];
            [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler_current, [0 0 0], [0 0 0], stressTensor, 'Mg', 'twin');
            twin_SFs = abs_schmid_factor(19:24,2);
            
            inds = ID==ID_current;
            
            A_grain = sum(inds(:));
            for iTwin = 1:6
                A_ref = sum(twin_map_ref(inds) == iTwin);
                A_def = sum(twin_map_def(inds) == iTwin);
                A_diff = A_ref - A_def;
                
                if A_ref > 0
                    vActiveTF = 1;
                else
                    vActiveTF = 0;
                end
                
                T3 = [T3; {iE_ref, iE, ID_current, iTwin, twin_SFs(iTwin), vActiveTF, A_grain, A_ref, A_def, A_diff}];
            end
        end
        
        ind =  (T3.vActiveTF==1)&(T3.iE==iE)&(T3.iE_ref==iE_ref); 
        t3 = T3(ind,:);
        edges = 0:0.05:0.5;
        vx = t3.variant_SF;
        vy = t3.A_diff ./ t3.A_grain;
        
        gv = discretize(vx, edges);
        nGroups = length(edges)-1;
        clear labels;
        for ii = 1:length(edges)-1
            labels{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
        end
        figure;disableDefaultInteractivity(gca);
        boxplot([vy; nan*ones(nGroups, 1)], [gv; (1:nGroups)'],'notch','on');   % prevent having empty group
        set(gca,'ylim',[-0.1 0.7]);
        
        xlabel('Twin Variant Schmid Factor');
        ylabel('Variant Detwin Area Fraction in Grain');
        set(gca,'xticklabels',labels,'xticklabelrotation',45,'fontsize',14);
        title(['load step ',num2str(iE), ' vs ',num2str(iE_ref), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','norma');
        print(fullfile(output_dir,['detwin fraction in grain vs SF iE_',num2str(iE),'_vs_',num2str(iE_ref),'.tiff']),'-dtiff');
        
        % The following uses variant fraction detwinned.
        vy = t3.A_diff ./ t3.A_ref;
        gv = discretize(vx, edges);
        nGroups = length(edges)-1;
        clear labels;
        for ii = 1:length(edges)-1
            labels{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
        end
        figure;disableDefaultInteractivity(gca);
        boxplot([vy; nan*ones(nGroups, 1)], [gv; (1:nGroups)'],'notch','on');   % prevent having empty group
        set(gca,'ylim',[-0.1 1.1]);
        
        xlabel('Twin Variant Schmid Factor');
        ylabel('Variant Detwin Fraction ');
        set(gca,'xticklabels',labels,'xticklabelrotation',45,'fontsize',14);
        title(['load step ',num2str(iE), ' vs ',num2str(iE_ref), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','norma');
        print(fullfile(output_dir,['detwin fraction vs SF iE_',num2str(iE),'_vs_',num2str(iE_ref),'.tiff']),'-dtiff');
        
        close all;
    end
    
end
%% Analyze 'significantly' detwin
% 'significantly detwinned variant': variant has completely detwinned 'grain',
% and accumulatively these completely detwinned grains accounts for >90% of
% the twinned area.
% (f1) variants partially/completely detwinned vs. SF
% also illustrate the 'past-or-present grains' belonging to the
% significantly detwinned variants

variableNames4 = {'iE_ref','iE','ID','iTwin','variant_SF', 'A_twin','A_residual','A_detwin', 'A_ratio'};
T4 = cell2table(cell(0,length(variableNames4)));
T4.Properties.VariableNames = variableNames4;

for iE_ref = 3
    iB_ef = iE_ref + 1;
    variant_ref = variant_cell{iB_ref};
    
    for iE = 4:6
        iB = iE + 1;
        variant_residual = variant_cell{iB};
        
        for ii = 1:length(gList)
            ID_current = gList(ii);
            ind_g = find(ID_current==gID);
            
            euler_current = [gPhi1(ind_g),gPhi(ind_g),gPhi2(ind_g)];
            [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler_current, [0 0 0], [0 0 0], stressTensor, 'Mg', 'twin');
            twin_SFs = abs_schmid_factor(19:24,2);
            
            inds = ismember(ID, ID_current);
            
            % find child twin grains (from past_or_present_twin_grain_ID)
            children = unique(past_or_present_twin_grain_ID_cell{iB}(inds));
            children(children==0) = [];
            children(isnan(children)) = [];
            
            % determine if a variant is 'significantly detwinned':
            % find active variants at iE=3
            active_variants = unique(variant_ref(inds));
            active_variants(active_variants==0) = [];
            
            % find twin grain for each active variant
            % check if this variant is significantly detwinned. i.e., if this
            % variant area at iE=3 contains mainly completely detwinned grain label
            % at iE=6
            for iTwin = 1:6
                if ismember(iTwin, active_variants)
                    inds_variant_ref = ismember(ID, ID_current) & ismember(variant_ref, iTwin);
                    size_variant_ref = sum(inds_variant_ref(:));
                    
                    inds_variant_residual = ismember(ID, ID_current) & ismember(variant_residual, iTwin);
                    size_variant_residual = sum(inds_variant_residual(:));
                    
                    labels = past_or_present_twin_grain_label_cell{iB}(inds_variant_ref); % g5 = completely detwin grain
                    g5 = 5;
                    size_detwin_grain = sum(labels==g5);
                    
                    if size_detwin_grain/size_variant_ref > 0.9
                        % significantly detwinned
                        T4 = [T4; {iE_ref, iE, ID_current, iTwin, twin_SFs(iTwin), size_variant_ref, size_variant_residual, size_detwin_grain, size_detwin_grain/size_variant_ref}];
                    else
                        T4 = [T4; {iE_ref, iE, ID_current, iTwin, twin_SFs(iTwin), size_variant_ref, size_variant_residual, size_detwin_grain, size_detwin_grain/size_variant_ref}];
                    end
                end
            end
            
        end
        
    end
    
end
%% [A1, Fig 15] Counts of significantly detwinned (and not) vs. twinSF 
for iE = 6 % 4:6
    iB = iE + 1;
    edges = -0.5:0.05:0.5;
    
    significantly_detwin_criterion = 0.9;
    ind1 = (T4.iE==iE) & (T4.A_ratio < significantly_detwin_criterion);    % 1: not significantly detwinned
    ind2 = (T4.iE==iE) & (T4.A_ratio >= significantly_detwin_criterion);
    
    N1 = histcounts(T4.variant_SF(ind1), edges);
    N2 = histcounts(T4.variant_SF(ind2), edges);
    
    figure; hold on; disableDefaultInteractivity(gca);
    hbar = bar(edges(1:end-1)+0.025, [N1(:), N2(:)], 1, 'stacked');
    
    hbar(1).FaceColor = [0 0 1];
    hbar(2).FaceColor = [1 0 0];
    
    set(gca,'fontsize',12, 'XTick',-0.5:0.1:0.5,'ylim',[0 150], 'xlim',[0 0.5]);
    xlabel('Twin Variant Schmid Factor');
    ylabel('Counts');
    
    yyaxis right;
    set(gca, 'ycolor', 'k','fontsize',16);
    ylabel('Percent (%)');
    plot(-0.475:0.05:0.475, N2./(N1+N2) * 100,'-ko','linewidth',1.5);
    set(gca,'ylim',[0 100]);
    
    % title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','norma');
    title(' ');
    legend({'Variants not significantly detwinned', 'Variants significantly detwinned','Percent of significantly detwinned'},'Location','northwest');
    print(fullfile(output_dir,['Fig 15 significantly detwin variant summary iE=',num2str(iE),'.tiff']), '-dtiff');
    
    
    % plot to check. label the grains with 'significantly detwinned variants'
    map = zeros(nR,nC);
    inds = past_or_present_twin_grain_label_cell{iB}>0;
    map(inds) = past_or_present_twin_grain_label_cell{iB}(inds);
    inds = current_twin_grain_label_cell{iB}>0;
    map(inds) = current_twin_grain_label_cell{iB}(inds);
    
    [f,a,c] = myplot(X,Y, map,  boundaryTF);
    
    colors = plasma(25);
    cmap = zeros(6,3);
    cmap(1,:) = [1,1,1];    % 0 = untwinned, background
    cmap(2,:) = colors(24,:);    % g1 = completely new twin, yellowish
    cmap(3,:) = colors(10,:);     % g2 = de-twin then re-twin, purple
    cmap(4,:) = colors(16,:);     % g3 = evoling current twin, pink
    cmap(5,:) = [0.7, 0.7, 0.7];     % g4 = evolving past-or-present twin, light gray
    cmap(6,:) = [0.2, 0.2, 0.2];    % g5 = completely de-twin, dark gray
    
    colormap(cmap);
    caxis([-0.5, 5.5]);
    set(c,'limits',[0.5, 5.5], 'Ticks', 1:5, ...
        'TickLabels', {'g1:fresh','g2:de-twin then re-twin','g3:evolving current','g4:evolving past or present','g5:completely detwin'})
    title(['grain summary, load step = ',num2str(iE), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','norma');
    title(' ');
    set(gcf,'position',[80,172,1000,600]);
    
    ids = T4.ID(ind2);
    label_map_with_ID(X,Y,ID,gcf,ids,'b',12,4);
    print(fullfile(output_dir,['Fig 15 significantly detwin variant grain iE=',num2str(iE),'.tiff']), '-dtiff');
    
    close all;
end


%% [Fig 3] plot strain map, twin variant map
for iE = 0:6
    iB = iE + 1;
    myplot(X,Y,variant_cell{iB},boundaryTF);
    set(gca,'xTickLabel',[],'yTickLabel',[]);
    title(['variant map, load step = ',num2str(iE), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','norma');
    print(fullfile(output_dir, ['variant_map_iE=',num2str(iE),'.tiff']), '-dtiff');
    close;
    
    d = matfile(fullfile(dic_dir, ['_',num2str(iE),'_v73.mat']));
    exx = d.exx;
    exx = exx(1:reduce_ratio:end,1:reduce_ratio:end);
    [f,a,c] = myplot(X,Y,exx,boundaryTF);
    caxis([-0.08 0.01]);
    set(gca,'xTickLabel',[],'yTickLabel',[], 'fontsize',16);
    % title(['exx map, load step = ',num2str(iE), ', \epsilon = ',num2str(strains(iB),'%.4f')],'fontweight','normal');
    title(['\epsilon^G = ',num2str(strains(iB),'%.4f')],'fontweight','normal');
    title(c,'\epsilon_{xx}')
    print(fullfile(output_dir, ['Fig 3 exx_map_iE=',num2str(iE),'.tiff']), '-dtiff');
    close;
end

%% plot basal slip traces for iE=1


%%
save(fullfile(output_dir,'temp_WS.mat'),'-v7.3');
% load(fullfile(output_dir,'temp_WS.mat'),'-v7.3');
