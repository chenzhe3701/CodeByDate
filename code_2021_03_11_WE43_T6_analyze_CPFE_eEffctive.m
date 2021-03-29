addChenFunction;
clear;
close all;
clc;
working_dir = 'E:\zhec umich Drive\WE43_T6_C1 EffStrain';

%% combine eff from CPFE
% d1 = load(fullfile(working_dir, 'EffStrain_Part2.mat'));    % for iE=1
% d = load(fullfile(working_dir, 'EffStrain.mat'));   % for iEs 2,3,4,5
% EffStrain_fromCPFE_Global = [d1.EffStrain_fromCPFE_Global(:); d.EffStrain_fromCPFE_Global(:)];
% X_new = d.X_new;
% Y_new = d.Y_new;
% ID_fromCPFE = d.ID_fromCPFE;
% save(fullfile(working_dir, 'EffStrainCombined_iE_1to5.mat'), 'EffStrain_fromCPFE_Global','ID_fromCPFE','X_new','Y_new');
% eEff = [d1.EffStrain_fromCPFE_Global(:)', d.EffStrain_fromCPFE_Global(:)'];

%% (1) load data from CPFE data
d = load(fullfile(working_dir, 'EffStrainCombined_iE_1to5.mat'));
eEff = d.EffStrain_fromCPFE_Global;
ID_CPFE = d.ID_fromCPFE;
X_CPFE = d.X_new;
Y_CPFE = d.Y_new;
boundary_CPFE = find_one_boundary_from_ID_matrix(ID_CPFE);
% myplot(X_CPFE, Y_CPFE, ID_CPFE, boundary_CPFE);

% find the same area from the adjusted map, to compare
saveDataPath = 'D:\WE43_T6_C1\Analysis_by_Matlab_after_realign';
load(fullfile(saveDataPath, 'WE43_T6_C1_EbsdToSemForTraceAnalysis_GbAdjusted'), ...
    'X','Y','boundaryTFB', 'uniqueBoundary','uniqueBoundaryList', ...
    'ID','gID','gNeighbors','gNNeighbors');

indR_min = find(Y(:,1) == Y_CPFE(1,1));
indR_max = find(Y(:,1) == Y_CPFE(end,1));
indC_min = find(X(1,:) == X_CPFE(1,1));
indC_max = find(X(1,:) == X_CPFE(1,end));
ind_step = (X_CPFE(1,2)-X_CPFE(1,1))/(X(1,2) - X(1,1));

X = X(indR_min:ind_step:indR_max, indC_min:ind_step:indC_max);
Y = Y(indR_min:ind_step:indR_max, indC_min:ind_step:indC_max);
ID = ID(indR_min:ind_step:indR_max, indC_min:ind_step:indC_max);
boundary = find_one_boundary_from_ID_matrix(ID);

% myplot(X, Y, ID, boundary);

d = load('D:\p\m\DIC_Analysis\temp_results\WE43_T6_C1_new_variant_map_20200401.mat','variantMapCleanedCell');
variant_map_cell = d.variantMapCleanedCell;

%% DIC strain data
dic_path = 'D:\WE43_T6_C1\SEM Data\stitched_DIC';

for iE = 1:5
    strainFile = matfile(fullfile(dic_path, ['_',num2str(iE),'_v73.mat']));
    eMap = calculate_effective_strain(strainFile.exx, strainFile.exy, strainFile.eyy);
    eMap_cell{iE} = eMap(indR_min:ind_step:indR_max, indC_min:ind_step:indC_max);   % crop
end

%% Now, get grains, and generating 9 types of distance maps, at each iE
x_dist = 1:250;

% Loop through iEs
for iE = 2:5
    % Load true twin map
    variant_map = variant_map_cell{iE}; % just for reference
    
    % find [bNum, gNum] pairs list, for different categories of activities
    % (1) not involved, (2) slip-induced-twin, slip side, (3) slip-induced-twin, twin side,
    % (4) co-found, (5) twin-induced-twin, old twin side, (6) twin-induced-twin, new twin side,
    % (7) slip-induced-growth, slip side, (8) slip-induced-growth, twin side, (9) co-growth
    
    % the data was stored in
    load(['D:\p\m\DIC_Analysis\temp_results\twin_gb_summary_',num2str(iE),'.mat'], ...
        'bg_not_involved', 'bg_slip_twin_a', 'bg_slip_twin_b', ...
        'bg_co_found','bg_twin_twin_a','bg_twin_twin_b', ...
        'bg_slip_growth_a','bg_slip_growth_b','bg_co_growth');
    
    eMap = eEff{iE-1};  % iE=1 is stored in cell{1}. but for iE=2, we need to look at eEff_iE_1, stored at cell{1}
    edmat = [];
    
    gID_list = unique(ID_CPFE(:));
    for ii = 1:length(gID_list)
        ID_current = gID_list(ii);
        ind = find(gID==ID_current);
        nNeighbors = gNNeighbors(ind);
        ID_neighbors = gNeighbors(ind, 1:nNeighbors);
        
        ind_local = ismember(ID_CPFE, [ID_current, ID_neighbors]);
        % Make it one data point wider on each side
        indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
        indC_max = min(size(ID_CPFE,2), find(sum(ind_local, 1), 1, 'last')+1);
        indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
        indR_max = min(size(ID_CPFE,1), find(sum(ind_local, 2), 1, 'last')+1);
        
        ID_local = ID_CPFE(indR_min:indR_max, indC_min:indC_max);
        eMap_local = eMap(indR_min:indR_max, indC_min:indC_max);
        
        for iNb = 1:nNeighbors
            ID_neighbor = ID_neighbors(iNb);
            
            % if neighbor grain in this local area
            if ismember(ID_neighbor, gID_list)
                % (1.1) Calculate this_uniqueGB number.
                if ID_current > ID_neighbor
                    gb = ID_current * 10000 + ID_neighbor;
                else
                    gb = ID_neighbor * 10000 + ID_current;
                end
                
                distMap_local = distance_from_boundary_in_grain(ID_local, [gb,ID_current]);
                % [[[1]]] Here is to summarize each grain boundary zone
                % individually, so that a point can be summarized multiple
                % times if it contribute to multiple gb zones
                edmat = [edmat; iE, ID_current, gb, arrayfun(@(x) nanmean(eMap_local(distMap_local==x)), 1:250)];   % recorded in x_dist=1:250
                
            end
        end
    end
    
    % We can calculate the distMap where the 'boundary grain list' can include all pairs    
    distMap_not_involved = distance_from_boundary_in_grain(ID_CPFE, bg_not_involved);
    distMap_slip_twin_a = distance_from_boundary_in_grain(ID_CPFE, bg_slip_twin_a);
    distMap_slip_twin_b = distance_from_boundary_in_grain(ID_CPFE, bg_slip_twin_b);
    distMap_co_found = distance_from_boundary_in_grain(ID_CPFE, bg_co_found);
    distMap_twin_twin_a = distance_from_boundary_in_grain(ID_CPFE, bg_twin_twin_a);
    distMap_twin_twin_b = distance_from_boundary_in_grain(ID_CPFE, bg_twin_twin_b);
    distMap_slip_growth_a = distance_from_boundary_in_grain(ID_CPFE, bg_slip_growth_a);
    distMap_slip_growth_b = distance_from_boundary_in_grain(ID_CPFE, bg_slip_growth_b);
    distMap_co_growth = distance_from_boundary_in_grain(ID_CPFE, bg_co_growth);
    
    % Plot, note: nanmean
    yLimits = [NaN       NaN;
        0.008    0.012;
        0.008    0.017;
        0.010    0.027;
        0.012    0.045];
    
    titleStr = {'','\epsilon^G = -0.004','\epsilon^G = -0.012','\epsilon^G = -0.023','\epsilon^G = -0.039'};
    
    um_per_dp = (X(1,2)-X(1,1)) * 360/4096;    % micron per data point, 88nm per pixel. DIC step size = 5 pxl, CPFE step size = 30 pixel. 
    nPts = 30;
    pt_range = 2:nPts;
    
    clear eline_not_involved eline_slip_twin_a eline_slip_twin_b eline_co_found eline_twin_twin_a eline_twin_twin_b eline_slip_growth_a eline_slip_growth_b eline_co_growth
    legend_str = []; istr = 1;
    
    figure;disableDefaultInteractivity(gca); hold on;
    
    % Note, edmat = [iE, ID, gb, es......]
    % This is the same method as in the paper. First calculate
    % f_eMean(dist_to_gb), then calculate mean of these data of grains
    % within each category.
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_not_involved, 'rows');
        eline_not_involved = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['not-involved: ',num2str(sum(ind))];
        legend_str{istr} = ['(1): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_not_involved(pt_range),'-','color',[0 0.6 0],'LineWidth',3);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_twin_a, 'rows');
        eline_slip_twin_a = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['slip-twin slip-side: ',num2str(sum(ind))];
        legend_str{istr} = ['(2): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_slip_twin_a(pt_range),'-b','LineWidth',3);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_growth_a, 'rows');
        eline_slip_growth_a = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['slip-growth slip-side: ',num2str(sum(ind))];
        legend_str{istr} = ['(3): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_slip_growth_a(pt_range),'-m','LineWidth',3);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_twin_b, 'rows');
        eline_slip_twin_b = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['slip-twin new-twin-side: ',num2str(sum(ind))];
        legend_str{istr} = ['(4): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_slip_twin_b(pt_range),'--b','LineWidth',2);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_co_found, 'rows');
        eline_co_found = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['co-found: ',num2str(sum(ind))];
        legend_str{istr} = ['(5): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_co_found(pt_range),'-k','LineWidth',2);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_twin_twin_b, 'rows');
        eline_twin_twin_b = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['twin-twin new-twin-side: ',num2str(sum(ind))];
        legend_str{istr} = ['(6): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_twin_twin_b(pt_range),'--r','LineWidth',2);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_growth_b, 'rows');
        eline_slip_growth_b = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['slip-growth twin-side: ',num2str(sum(ind))];
        legend_str{istr} = ['(7): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_slip_growth_b(pt_range),'--m','LineWidth',.5);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_twin_twin_a, 'rows');
        eline_twin_twin_a = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['twin-twin old-twin-side: ',num2str(sum(ind))];
        legend_str{istr} = ['(8): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_twin_twin_a(pt_range),'-r','LineWidth',.5);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_co_growth, 'rows');
        eline_co_growth = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['co-growth: ',num2str(sum(ind))];
        legend_str{istr} = ['(9): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_co_growth(pt_range),'-c','LineWidth',.5);
    end
    
    legend(legend_str, 'location','eastoutside');
    % set(gca,'fontsize',16,'ylim',yLimits(iE,:));
    xlabel(['Distance to Grain Boundary (',char(181),'m)']);
    ylabel('Mean of Effective Strain');
    
    title(['iE=',num2str(iE)],'fontweight','normal');
    title(titleStr{iE},'fontweight','normal');
    yLimits(iE,:) = get(gca,'YLim');
    print(fullfile(working_dir, ['cpfe_gAvg_iE_',num2str(iE),'.tiff']), '-dtiff');
    close;
    
    % [[[2]]] try another method, directly summarize from dist map
    if 1
        clear edmat_cat
        edmat_cat(1,:) = arrayfun(@(x) nanmean(eMap(distMap_not_involved==x)), 1:250);
        edmat_cat(2,:) = arrayfun(@(x) nanmean(eMap(distMap_slip_twin_a==x)), 1:250);
        edmat_cat(3,:) = arrayfun(@(x) nanmean(eMap(distMap_slip_growth_a==x)), 1:250);
        
        edmat_cat(4,:) = arrayfun(@(x) nanmean(eMap(distMap_slip_twin_b==x)), 1:250);
        edmat_cat(5,:) = arrayfun(@(x) nanmean(eMap(distMap_co_found==x)), 1:250);
        edmat_cat(6,:) = arrayfun(@(x) nanmean(eMap(distMap_twin_twin_b==x)), 1:250);
        
        edmat_cat(7,:) = arrayfun(@(x) nanmean(eMap(distMap_slip_growth_b==x)), 1:250);
        edmat_cat(8,:) = arrayfun(@(x) nanmean(eMap(distMap_twin_twin_a==x)), 1:250);
        edmat_cat(9,:) = arrayfun(@(x) nanmean(eMap(distMap_co_growth==x)), 1:250);
        
        figure; hold on;
        istr = 1;
        try
            legend_str{istr} = ['(1) not-involved: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(1,pt_range),'-','color',[0 0.6 0],'LineWidth',3);
        end
        
        try
            legend_str{istr} = ['(2) slip-twin slip-side: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(2,pt_range),'-b','LineWidth',3);
        end
        
        try
            legend_str{istr} = ['(3) slip-growth slip-side: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(3,pt_range),'-m','LineWidth',3);
        end
        
        try
            legend_str{istr} = ['(4) slip-twin new-twin-side: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(4,pt_range),'--b','LineWidth',2);
        end
        
        try
            legend_str{istr} = ['(5) co-found: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(5,pt_range),'-k','LineWidth',2);
        end
        
        try
            legend_str{istr} = ['(6) twin-twin new-twin-side: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(6,pt_range),'--r','LineWidth',2);
        end
        
        try
            legend_str{istr} = ['(7) slip-growth twin-side: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(7,pt_range),'--m','LineWidth',.5);
        end
        
        try
            legend_str{istr} = ['(8) twin-twin old-twin-side: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(8,pt_range),'-r','LineWidth',.5);
        end
        
        try
            legend_str{istr} = ['(9) co-growth: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(9,pt_range),'-c','LineWidth',.5);
        end
        
        legend(legend_str, 'location','eastoutside');
        % set(gca,'fontsize',16,'ylim',yLimits(iE,:));
        xlabel(['Distance to Grain Boundary (',char(181),'m)']);
        ylabel('Mean of Effective Strain');
        
        title(['iE=',num2str(iE)],'fontweight','normal');
        title(titleStr{iE},'fontweight','normal');
        yLimits(iE,:) = get(gca,'YLim');
        
        print(fullfile(working_dir, ['cpfe_pAvg_iE_',num2str(iE),'.tiff']), '-dtiff');
        close;
    end
    
    save(fullfile(working_dir,['data_iE=',num2str(iE),'.mat']), 'edmat', 'edmat_cat', ...
        'distMap_not_involved', 'distMap_slip_twin_a', 'distMap_slip_twin_b', ...
        'distMap_co_found', 'distMap_twin_twin_a', 'distMap_twin_twin_b', ...
        'distMap_slip_growth_a', 'distMap_slip_growth_b', 'distMap_co_growth');
    
    
end

%% [[for comparison]] Using experimental strain of the same CPFE area

x_dist = 1:250;

% Loop through iEs
for iE = 2:5
    % Load true twin map
    variant_map = variant_map_cell{iE}; % just for reference
    
    % find [bNum, gNum] pairs list, for different categories of activities
    % (1) not involved, (2) slip-induced-twin, slip side, (3) slip-induced-twin, twin side,
    % (4) co-found, (5) twin-induced-twin, old twin side, (6) twin-induced-twin, new twin side,
    % (7) slip-induced-growth, slip side, (8) slip-induced-growth, twin side, (9) co-growth
    
    % the data was stored in
    load(['D:\p\m\DIC_Analysis\temp_results\twin_gb_summary_',num2str(iE),'.mat'], ...
        'bg_not_involved', 'bg_slip_twin_a', 'bg_slip_twin_b', ...
        'bg_co_found','bg_twin_twin_a','bg_twin_twin_b', ...
        'bg_slip_growth_a','bg_slip_growth_b','bg_co_growth');
        
    eMap = eMap_cell{iE-1};     % =====================> This is the important part. eMap is from exprimental results, rather than CPFE.   

    edmat = [];
    
    gID_list = unique(ID(:));
    for ii = 1:length(gID_list)
        ID_current = gID_list(ii);
        ind = find(gID==ID_current);
        nNeighbors = gNNeighbors(ind);
        ID_neighbors = gNeighbors(ind, 1:nNeighbors);
        
        ind_local = ismember(ID, [ID_current, ID_neighbors]);
        % Make it one data point wider on each side
        indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
        indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
        indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
        indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        eMap_local = eMap(indR_min:indR_max, indC_min:indC_max);
        
        for iNb = 1:nNeighbors
            ID_neighbor = ID_neighbors(iNb);
            
            % if neighbor grain in this local area
            if ismember(ID_neighbor, gID_list)
                % (1.1) Calculate this_uniqueGB number.
                if ID_current > ID_neighbor
                    gb = ID_current * 10000 + ID_neighbor;
                else
                    gb = ID_neighbor * 10000 + ID_current;
                end
                
                distMap_local = distance_from_boundary_in_grain(ID_local, [gb,ID_current]);
                edmat = [edmat; iE, ID_current, gb, arrayfun(@(x) nanmean(eMap_local(distMap_local==x)), 1:250)];   % recorded in x_dist=1:250
                
            end
        end
    end
    
    % We can calculate the distMap where the 'boundary grain list' can include all pairs
    distMap_not_involved = distance_from_boundary_in_grain(ID_CPFE, bg_not_involved);
    distMap_slip_twin_a = distance_from_boundary_in_grain(ID_CPFE, bg_slip_twin_a);
    distMap_slip_twin_b = distance_from_boundary_in_grain(ID_CPFE, bg_slip_twin_b);
    distMap_co_found = distance_from_boundary_in_grain(ID_CPFE, bg_co_found);
    distMap_twin_twin_a = distance_from_boundary_in_grain(ID_CPFE, bg_twin_twin_a);
    distMap_twin_twin_b = distance_from_boundary_in_grain(ID_CPFE, bg_twin_twin_b);
    distMap_slip_growth_a = distance_from_boundary_in_grain(ID_CPFE, bg_slip_growth_a);
    distMap_slip_growth_b = distance_from_boundary_in_grain(ID_CPFE, bg_slip_growth_b);
    distMap_co_growth = distance_from_boundary_in_grain(ID_CPFE, bg_co_growth);
    
    % Plot, note: nanmean
    yLimits = [NaN       NaN;
        0.008    0.012;
        0.008    0.017;
        0.010    0.027;
        0.012    0.045];
    
    titleStr = {'','\epsilon^G = -0.004','\epsilon^G = -0.012','\epsilon^G = -0.023','\epsilon^G = -0.039'};
    
    um_per_dp = (X(1,2)-X(1,1)) * 360/4096;    % micron per data point, ~0.43
    
    nPts = 30;
    pt_range = 2:nPts;
    
    clear eline_not_involved eline_slip_twin_a eline_slip_twin_b eline_co_found eline_twin_twin_a eline_twin_twin_b eline_slip_growth_a eline_slip_growth_b eline_co_growth
    legend_str = []; istr = 1;
    
    figure;disableDefaultInteractivity(gca); hold on;
    
    % Note, edmat = [iE, ID, gb, es......]
    % This is the same method as in the paper. First calculate
    % f_eMean(dist_to_gb), then calculate mean of these data of grains
    % within each category.
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_not_involved, 'rows');
        eline_not_involved = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['not-involved: ',num2str(sum(ind))];
        legend_str{istr} = ['(1): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_not_involved(pt_range),'-','color',[0 0.6 0],'LineWidth',3);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_twin_a, 'rows');
        eline_slip_twin_a = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['slip-twin slip-side: ',num2str(sum(ind))];
        legend_str{istr} = ['(2): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_slip_twin_a(pt_range),'-b','LineWidth',3);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_growth_a, 'rows');
        eline_slip_growth_a = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['slip-growth slip-side: ',num2str(sum(ind))];
        legend_str{istr} = ['(3): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_slip_growth_a(pt_range),'-m','LineWidth',3);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_twin_b, 'rows');
        eline_slip_twin_b = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['slip-twin new-twin-side: ',num2str(sum(ind))];
        legend_str{istr} = ['(4): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_slip_twin_b(pt_range),'--b','LineWidth',2);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_co_found, 'rows');
        eline_co_found = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['co-found: ',num2str(sum(ind))];
        legend_str{istr} = ['(5): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_co_found(pt_range),'-k','LineWidth',2);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_twin_twin_b, 'rows');
        eline_twin_twin_b = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['twin-twin new-twin-side: ',num2str(sum(ind))];
        legend_str{istr} = ['(6): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_twin_twin_b(pt_range),'--r','LineWidth',2);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_growth_b, 'rows');
        eline_slip_growth_b = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['slip-growth twin-side: ',num2str(sum(ind))];
        legend_str{istr} = ['(7): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_slip_growth_b(pt_range),'--m','LineWidth',.5);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_twin_twin_a, 'rows');
        eline_twin_twin_a = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['twin-twin old-twin-side: ',num2str(sum(ind))];
        legend_str{istr} = ['(8): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_twin_twin_a(pt_range),'-r','LineWidth',.5);
    end
    
    try
        [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_co_growth, 'rows');
        eline_co_growth = nanmean(edmat(ind,4:end),1);
        legend_str{istr} = ['co-growth: ',num2str(sum(ind))];
        legend_str{istr} = ['(9): ',num2str(sum(ind))];
        istr = istr + 1;
        plot(um_per_dp*x_dist(pt_range),eline_co_growth(pt_range),'-c','LineWidth',.5);
    end
    
    legend(legend_str, 'location','eastoutside');
    % set(gca,'fontsize',16,'ylim',yLimits(iE,:));
    xlabel(['Distance to Grain Boundary (',char(181),'m)']);
    ylabel('Mean of Effective Strain');
    
    title(['iE=',num2str(iE)],'fontweight','normal');
    title(titleStr{iE},'fontweight','normal');
    yLimits(iE,:) = get(gca,'YLim');
    
    print(fullfile(working_dir, ['exp_gAvg_iE_',num2str(iE),'.tiff']), '-dtiff');
    close;
    
    % [[[2]]] try another method, directly summarize from dist map
    if 1
        clear edmat_cat
        edmat_cat(1,:) = arrayfun(@(x) nanmean(eMap(distMap_not_involved==x)), 1:250);
        edmat_cat(2,:) = arrayfun(@(x) nanmean(eMap(distMap_slip_twin_a==x)), 1:250);
        edmat_cat(3,:) = arrayfun(@(x) nanmean(eMap(distMap_slip_growth_a==x)), 1:250);
        
        edmat_cat(4,:) = arrayfun(@(x) nanmean(eMap(distMap_slip_twin_b==x)), 1:250);
        edmat_cat(5,:) = arrayfun(@(x) nanmean(eMap(distMap_co_found==x)), 1:250);
        edmat_cat(6,:) = arrayfun(@(x) nanmean(eMap(distMap_twin_twin_b==x)), 1:250);
        
        edmat_cat(7,:) = arrayfun(@(x) nanmean(eMap(distMap_slip_growth_b==x)), 1:250);
        edmat_cat(8,:) = arrayfun(@(x) nanmean(eMap(distMap_twin_twin_a==x)), 1:250);
        edmat_cat(9,:) = arrayfun(@(x) nanmean(eMap(distMap_co_growth==x)), 1:250);
        
        figure; hold on;
        istr = 1;
        try
            legend_str{istr} = ['(1) not-involved: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(1,pt_range),'-','color',[0 0.6 0],'LineWidth',3);
        end
        
        try
            legend_str{istr} = ['(2) slip-twin slip-side: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(2,pt_range),'-b','LineWidth',3);
        end
        
        try
            legend_str{istr} = ['(3) slip-growth slip-side: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(3,pt_range),'-m','LineWidth',3);
        end
        
        try
            legend_str{istr} = ['(4) slip-twin new-twin-side: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(4,pt_range),'--b','LineWidth',2);
        end
        
        try
            legend_str{istr} = ['(5) co-found: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(5,pt_range),'-k','LineWidth',2);
        end
        
        try
            legend_str{istr} = ['(6) twin-twin new-twin-side: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(6,pt_range),'--r','LineWidth',2);
        end
        
        try
            legend_str{istr} = ['(7) slip-growth twin-side: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(7,pt_range),'--m','LineWidth',.5);
        end
        
        try
            legend_str{istr} = ['(8) twin-twin old-twin-side: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(8,pt_range),'-r','LineWidth',.5);
        end
        
        try
            legend_str{istr} = ['(9) co-growth: '];
            istr = istr + 1;
            plot(um_per_dp*x_dist(pt_range),edmat_cat(9,pt_range),'-c','LineWidth',.5);
        end
        
        legend(legend_str, 'location','eastoutside');
        % set(gca,'fontsize',16,'ylim',yLimits(iE,:));
        xlabel(['Distance to Grain Boundary (',char(181),'m)']);
        ylabel('Mean of Effective Strain');
        
        title(['iE=',num2str(iE)],'fontweight','normal');
        title(titleStr{iE},'fontweight','normal');
        yLimits(iE,:) = get(gca,'YLim');
        
        print(fullfile(working_dir, ['exp_pAvg_iE_',num2str(iE),'.tiff']), '-dtiff');
        close;
    end
    
    save(fullfile(working_dir,['exp_data_iE=',num2str(iE),'.mat']), 'edmat', 'edmat_cat', ...
        'distMap_not_involved', 'distMap_slip_twin_a', 'distMap_slip_twin_b', ...
        'distMap_co_found', 'distMap_twin_twin_a', 'distMap_twin_twin_b', ...
        'distMap_slip_growth_a', 'distMap_slip_growth_b', 'distMap_co_growth');
    
    
end










