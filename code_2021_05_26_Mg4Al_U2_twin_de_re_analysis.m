
%% Summarize the twin-detwin-retwin behavior for Mg4Al_U2 (should be U2)
clear; clc; close all;
addChenFunction;

sample_name = 'Mg4Al_U2';

% location of the parent/twin grain/variant/... files
data_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis';

% location of the twin-detwin-retwin data
data_dir_2 = 'E:\zhec umich Drive\0_temp_output\2021-05-04 Analyze persistent twin\Mg4Al_U2 reMake';
data = matfile(fullfile(data_dir_2, 'Mg4Al_U2 twin data.mat'));

% Where to save the output figures
output_dir = ['E:\zhec umich Drive\0_temp_output\2021-05-04 Analyze persistent twin\',sample_name, ' analysis reMake'];
mkdir(output_dir);

strain_ebsd = [0, -0.0085, -0.0220, -0.0342, ...
    -0.0343, -0.0268, -0.0150, -0.0032, ...
    -0.0112, -0.0207, -0.0297, ...
    -0.0235, -0.0130, -0.0045];
%% load reference data at iE = 0, reference ID map and grain boundary
iE = 0;
d = matfile(fullfile(data_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
ID_0 = d.ID;
gPhi1_0 = d.gPhi1;
gPhi_0 = d.gPhi;
gPhi2_0 = d.gPhi2;
gDiameter_0 = d.gDiameter;
gID_0 = d.gID;
x = d.x;
y = d.y;

boundary_0 = find_one_boundary_from_ID_matrix(ID_0);

%% load from data file
ID_p_iB_to_1_cell = data.ID_p_iB_to_1_cell;
ID_c_iB_to_1_cell = data.ID_c_iB_to_1_cell;
variant_pixel_iB_to_1_cell = data.variant_pixel_iB_to_1_cell;

past_or_present_twin_cell = data.past_or_present_twin_cell;     
current_twin_cell = data.current_twin_cell;

past_or_present_twin_grain_ID_cell = data.past_or_present_twin_grain_ID_cell;
past_or_present_twin_grain_label_cell = data.past_or_present_twin_grain_label_cell;
current_twin_grain_label_cell = data.current_twin_grain_label_cell;

variant_pixel_cell = data.variant_pixel_cell;
boundary_p_iB_to_1_cell = data.boundary_p_iB_to_1_cell;
variant_pixel_mode = data.variant_pixel_mode;

gList = data.gList;
ID_overlap = data.ID_overlap;

%% Reconstruct variant map
% Generate a variant map for each iE, using information from all iEs:
% The variant# is the most frequent (the 'mode') variant at this pixel
% position among all iEs.
% the 'past_or_present map might contain more twin pixels, so we need to
% find the 'valid area' is the twin area at this iE
for iE = 0:13
    iB = iE + 1;
    
    ind_past_or_present = past_or_present_twin_cell{iB} > 0; 
    map = variant_pixel_mode;
    map(~ind_past_or_present) = 0;
    combined_past_or_present_variant_cell{iB} = map;
    
    ind_current_twin = current_twin_cell{iB}>0;
    map = variant_pixel_mode;
    map(~ind_current_twin) = 0;
    combined_current_variant_cell{iB} = map;
    
    [f,a,c] = myplot(x, y, combined_past_or_present_variant_cell{iB}, boundary_p_iB_to_1_cell{iB});
    caxis([0,6]);
    title(['combined past or present variant, load step = ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    set(a, 'xTick',[], 'yTick',[]);
    print(fullfile(output_dir,['combined past or present variant iE=',num2str(iE),'.tiff']),'-dtiff');
    
    [f,a,c] = myplot(x, y, combined_current_variant_cell{iB}, boundary_p_iB_to_1_cell{iB});
    caxis([0,6]);
    title(['combined current variant, load step = ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    set(a, 'xTick',[], 'yTick',[]);
    print(fullfile(output_dir,['combined current variant iE=',num2str(iE),'.tiff']),'-dtiff');
    
    close all;
end



%% Summaries: Counts of twin/non-twin variants/grains vs. SF

variableNames = {'iE','ID','max_basal_SF','max_twin_SF','activeTwin_SF','gDia','twinnedTF','tAF'};
T = cell2table(cell(0,length(variableNames)));
T.Properties.VariableNames = variableNames;

variableNames2 = {'iE','ID','iTwin','variant','variant_SF','vActiveTF','vPct','max_basal_SF','max_twin_SF'};
T2 = cell2table(cell(0,length(variableNames2)));
T2.Properties.VariableNames = variableNames2;

stressTensor = [-1 0 0; 0 0 0; 0 0 0];

% For each parent grian, determine child grain
for iE = 0:13
    iB = iE + 1;  
    
    ID_p = ID_p_iB_to_1_cell{iB};
    ID_c = ID_c_iB_to_1_cell{iB};
    
    for ii = 1:length(gList)
        ID_p_current =  gList(ii);
        
        inds = ismember(ID_overlap, ID_p_current);  % ========> try to use ID_overlap  
        
        % Find area fraction for each variant
        vol_grain = sum(inds(:));
        variant_ids = combined_current_variant_cell{iB}(inds);
        vol_variants = arrayfun(@(x) sum(ismember(variant_ids, x)), 1:6);
        tAF = sum(vol_variants)/sum(vol_grain);
        
        % If variant is twinned at this strain level.
        activeTS = vol_variants>0;
        
        ind_g = find(ID_p_current==gID_0);
        gd = gDiameter_0(ind_g);
        
        euler_current = [gPhi1_0(ind_g),gPhi_0(ind_g),gPhi2_0(ind_g)];
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
        
        T = [T; {iE, ID_p_current, max_basal_SF, max_twin_SF, activeTwin_SF, gd, twinnedTF, tAF}];
        
        for iTwin=1:length(activeTS)
            if activeTS(iTwin)>0
                vActiveTF = true;
            else
                vActiveTF = false;
            end
            vi = SF_order(iTwin);
            vSF = twin_SFs(iTwin);
            vPct = vol_variants(iTwin)/vol_grain;
            T2 = [T2; {iE, ID_p_current, iTwin, vi, vSF, vActiveTF, vPct, max_basal_SF, max_twin_SF}];
        end
        
        if rem(ii,10) == 1
            str = sprintf('iE=%d, ID=%d', iE, ID_p_current);
            disp(str);
        end
        
        
    end
end


%% Fig 2a. Counts of variants twinned & not-twinned vs. variant_SF
edges = -0.5:0.05:0.5;
for iE = 0:13
    iB = iE + 1;
    ind = (T2.iE==iE)&(T2.vActiveTF == 1);
    ind2 = (T2.iE==iE)&(T2.vActiveTF == 0);
    [N_t,~] = histcounts(T2.variant_SF(ind), edges);
    [N_nt,~] = histcounts(T2.variant_SF(ind2), edges);
    
    figure; hold on; disableDefaultInteractivity(gca);
    hbar = bar(edges(1:end-1)+0.025, [N_nt(:), N_t(:)], 1, 'stacked');
    set(gca,'fontsize',12, 'XTick',-0.5:0.1:0.5,'ylim',[0 300]);
    xlabel('Twin Variant Schmid Factor');
    ylabel('Counts');
    
    yyaxis right;
    set(gca, 'ycolor', 'k','fontsize',16);
    ylabel('Percent (%)');
    plot(-0.475:0.05:0.475, N_t./(N_t+N_nt) * 100,'-ko','linewidth',1.5);
    ylim = get(gca,'ylim');
    set(gca, 'ylim',[0 100]);
    
    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    %     title(' ');
    legend({'Variants not twinned', 'Variants twinned','Percent of variants twinned'},'Location','northwest');
    
    hbar(1).FaceColor = [0 0 1];
    hbar(2).FaceColor = [1 0 0];
           
    print(fullfile(output_dir,['twin nontwin variant vs SF iE=',num2str(iE),'.tiff']),'-dtiff');

    disp('row1: # twinned, row2: # not twinned, row3: pct');
    clear tt;
    tt = [N_t; N_nt; N_t./(N_t+N_nt)];
    disp(array2table(tt));
    
    close all;
end


%% Fig 2b. Counts grains twinned & not twinned vs. max_basal_SF
edges = 0:0.05:0.5;
for iE = 0:13
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
    set(gca,'ylim',[0 100]);
    xlabel('Maximum Basal Schmid Factor');
    ylabel('Counts');
    
    yyaxis right;
    set(gca,'ycolor','k','XTick',edges(1:end-1)+d_int,'xTickLabels',xstr,'xTickLabelRotation',45);
    plot(xpos, N_t./(N_t+N_nt) * 100,'-ko','linewidth',1.5);
    ylabel('Percent (%)');
    legend('Grains not twinned','Grains twinned','Percent of grains twinned','location','northwest');
    set(gca,'fontsize',12,'ylim',[0 149],'fontsize',16);
    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    % title('');
    hbar(1).FaceColor = [0 0 1];
    hbar(2).FaceColor = [1 0 0];
        
    print(fullfile(output_dir,['twin nontwin grain vs basal SF iE=',num2str(iE),'.tiff']),'-dtiff');
    
    disp(['iE = ',num2str(iE)]);
    tt = [N_t; N_nt; N_t./(N_t+N_nt)];
    disp(array2table(tt'));
    
    close all;
end


% Fig 2c/d
for iE = 0:13
    iB = iE + 1;
    figure;disableDefaultInteractivity(gca); hold on;
    ind = (T.iE==iE)&(T.twinnedTF==0);
    plot(T.max_twin_SF(ind), T.max_basal_SF(ind),'.','markersize',12,'color','b'); % '#0072BD'
    xlabel('Max Twin Schmid Factor');
    ylabel('Max Basal Schmid Factor');
    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    legend('Not twinned','location','southwest');
    set(gca,'fontsize',16,'xTick',-0.5:0.1:0.5,'xlim',[-0.5,0.5],'ylim',[0 0.5]);
    % title('')

    figure;disableDefaultInteractivity(gca); hold on;
    ind = (T.iE==iE)&(T.twinnedTF==1);
    plot(T.max_twin_SF(ind), T.max_basal_SF(ind),'.','markersize',12,'color','r');    % '#D95319'
    xlabel('Max Twin Schmid Factor');
    ylabel('Max Basal Schmid Factor');
    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    legend('Twinned','location','southwest');
    set(gca,'fontsize',16,'xTick',-0.5:0.25:0.5,'xlim',[-0.5,0.5],'ylim',[0 0.5]);
    % title('');
    
    figure;disableDefaultInteractivity(gca); hold on;
    ind = (T.iE==iE)&(T.twinnedTF==1);
    plot(T.max_twin_SF(ind), T.max_basal_SF(ind),'.','markersize',12,'color','r');    % '#D95319'
    ind = (T.iE==iE)&(T.twinnedTF==0);
    plot(T.max_twin_SF(ind), T.max_basal_SF(ind),'.','markersize',12,'color','b'); % '#0072BD'
    xlabel('Max Twin Schmid Factor');
    ylabel('Max Basal Schmid Factor');
    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    legend('Twinned', 'Not twinned', 'location','southwest');
    set(gca,'fontsize',16,'xTick',-0.5:0.25:0.5,'xlim',[-0.5,0.5],'ylim',[0 0.5]);
    % title('');
    
    print(fullfile(output_dir,['twin nontwin grain vs max SF iE=',num2str(iE),'.tiff']),'-dtiff');
    
    close all;
end

%% Fig 3. variant_area_fraction_in_grain vs. variant_SF
for iE = 0:13
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
    boxplot([vy; nan*ones(nGroups, 1)], [gv; (1:nGroups)'],'notch','on');   % prevent having empty group
    
    xlabel('Twin Variant Schmid Factor');
    ylabel('Twin Variant Area Fraction in Grain');
    set(gca,'xticklabels',labels,'xticklabelrotation',45,'fontsize',15);
    title(['twinned grains, load step = ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    % title([' ']);
    print(fullfile(output_dir,['variant fraction vs SF iE=',num2str(iE),'.tiff']),'-dtiff');
    
    close all;
end

%% (A2) variant detwin area (as pct of twin area OR as pct of grain area?) vs. SF 
% This can use pixel map, just compare iE_ref (which has max twin area, such as iE=3) vs iE (which has decreased twin area)

variableNames3 = {'iE_ref','iE','ID','iTwin','variant_SF','vActiveTF','A_grain','A_ref','A_def','A_diff'};
T3 = cell2table(cell(0,length(variableNames3)));
T3.Properties.VariableNames = variableNames3;

for iE_ref = 3
    iB_ref = iE_ref + 1;
    twin_map_ref = combined_past_or_present_variant_cell{iB_ref};
    
    for iE = 4:7
        iB = iE + 1;
        twin_map_def = combined_current_variant_cell{iB};
        
        for ii = 1:length(gList)
            ID_current = gList(ii);
            
            ind_g = find(ID_current==gID_0);
            gd = gDiameter_0(ind_g);
            
            euler_current = [gPhi1_0(ind_g), gPhi_0(ind_g), gPhi2_0(ind_g)];
            [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler_current, [0 0 0], [0 0 0], stressTensor, 'Mg', 'twin');
            twin_SFs = abs_schmid_factor(19:24,2);
            
            inds = ismember(ID_overlap, ID_current);
            
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
                
                % ==> If variant area > 0 at iE_ref, then variant active. 
                % Maybe also interesting to check those with only 'A_def' > 0   
                T3 = [T3; {iE_ref, iE, ID_current, iTwin, twin_SFs(iTwin), vActiveTF, A_grain, A_ref, A_def, A_diff}];
            end
        end
        
        ind = (T3.vActiveTF==1)&(T3.iE==iE)&(T3.iE_ref==iE_ref);
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
        % set(gca,'ylim',[-0.01 1]);
        xlabel('Twin Variant Schmid Factor');
        ylabel('Variant Detwin Area Fraction in Grain');
        set(gca,'xticklabels',labels,'xticklabelrotation',45,'fontsize',14);
        title(['load step ',num2str(iE), ' vs ',num2str(iE_ref), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
        % title([' '],'fontweight','normal');
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
        % set(gca,'ylim',[-0.01 2]);
        xlabel('Twin Variant Schmid Factor');
        ylabel('Variant Detwin Fraction ');
        set(gca,'xticklabels',labels,'xticklabelrotation',45,'fontsize',14);
        set(gca,'ylim', [-0.1, 1.1]);
        title(['load step ',num2str(iE), ' vs ',num2str(iE_ref), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
        % title([' '],'fontweight','normal');
        print(fullfile(output_dir,['detwin fraction vs SF iE_',num2str(iE),'_vs_',num2str(iE_ref),'.tiff']),'-dtiff');
        
        close all;
        
    end
    
end


%% Analyze detwin
% Method-1: 'significantly detwinned variant': variant has completely
% detwinned 'grain', and accumulatively these completely detwinned grains
% accounts for >90% of the twinned area.

variableNames4 = {'iE_ref','iE','ID','iTwin','variant_SF', 'A_twin','A_residual','A_detwin', 'A_ratio'};
T4 = cell2table(cell(0,length(variableNames4)));
T4.Properties.VariableNames = variableNames4;

for iE_ref = 3
    iB_ref = iE_ref + 1;
    variant_ref = combined_past_or_present_variant_cell{iB_ref};
    
    for iE = 4:7
        iB = iE + 1;
        variant_residual = combined_current_variant_cell{iB};
        
        for ii = 1:length(gList)
            ID_current = gList(ii);
            ind_g = find(ID_current==gID_0);
            
            euler_current = [gPhi1_0(ind_g), gPhi_0(ind_g), gPhi2_0(ind_g)];
            [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler_current, [0 0 0], [0 0 0], stressTensor, 'Mg', 'twin');
            twin_SFs = abs_schmid_factor(19:24,2);
            
            inds = ismember(ID_overlap, ID_current);
            
            % find child twin grains (from past_or_present_twin_grain_ID)
            children = unique(past_or_present_twin_grain_ID_cell{iB}(inds));
            children(children==0) = [];
            children(isnan(children)) = [];
            
            % Method-1: determine if a variant is 'significantly detwinned':
            % find active variants at iE_ref (=3)
            active_variants = unique(variant_ref(inds));
            active_variants(active_variants==0) = [];
            
            % find twin grain for each active variant
            % check if this variant is significantly detwinned. i.e., if this
            % variant area at iE=3 contains mainly completely detwinned grain label at iE=6 
            for iTwin = 1:6
                if ismember(iTwin, active_variants)
                    inds_variant_ref = ismember(ID_overlap, ID_current) & ismember(variant_ref, iTwin);
                    size_variant_ref = sum(inds_variant_ref(:));
                    
                    inds_variant_residual = ismember(ID_overlap, ID_current) & ismember(variant_residual, iTwin);
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


% [A1] Counts of significantly detwinned (and not) vs. twinSF 
for iE = 4:7
    iB = iE + 1;
    
    edges = -0.5:0.05:0.5;
    ind1 = (T4.iE==iE) & (T4.A_ratio < 0.9);    % 1: NOT significantly detwinned
    ind2 = (T4.iE==iE) & (T4.A_ratio >= 0.9);
    N1 = histcounts(T4.variant_SF(ind1), edges);
    N2 = histcounts(T4.variant_SF(ind2), edges);
    
    figure; hold on; disableDefaultInteractivity(gca);
    hbar = bar(edges(1:end-1)+0.025, [N1(:), N2(:)], 1, 'stacked');
    
    hbar(1).FaceColor = [0 0 1];
    hbar(2).FaceColor = [1 0 0];
    
    set(gca,'fontsize',12, 'XTick',-0.5:0.1:0.5);
    ylim = get(gca,'ylim');
    set(gca, 'ylim', [0 1.5 * ylim(2)]);
    xlabel('Twin Variant Schmid Factor');
    ylabel('Counts');
    
    yyaxis right;
    set(gca, 'ycolor', 'k','fontsize',16);
    ylabel('Percent (%)');
    plot(-0.475:0.05:0.475, N2./(N1+N2) * 100,'-ko','linewidth',1.5);
    ylim = get(gca,'ylim');
    set(gca, 'ylim', [0 1.5 * ylim(2)]);
    
    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    legend({'Variants not significantly detwinned', 'Variants significantly detwinned','Percent of significantly detwinned'},'Location','northwest');
    print(fullfile(output_dir,['method-1 detwin variant vs SF iE=',num2str(iE),'.tiff']), '-dtiff');
    
    
    %% plot to check. label the grains with 'significantly detwinned variants'
    map = zeros(size(ID_0));
    inds = past_or_present_twin_grain_label_cell{iB}>0;
    map(inds) = past_or_present_twin_grain_label_cell{iB}(inds);
    inds = current_twin_grain_label_cell{iB}>0;
    map(inds) = current_twin_grain_label_cell{iB}(inds);
    
    [f,a,c] = myplot(x,y, map,  boundary_p_iB_to_1_cell{iB});
    
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
    title(['grain summary, load step = ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    set(gcf,'position',[80,172,1000,600]);
    
    % label those grains
    ids = T4.ID(ind2);
    label_map_with_ID(x,y,ID_overlap, gcf,ids,'r',12,4);
    print(fullfile(output_dir,['method-1 detwin variant grain iE',num2str(iE),'.tiff']), '-dtiff');
    
    close all;
end

%% Method-2:
% Contains 'completely detwinned grain', then the major variant of that grain is recorded

variableNames5 = {'iE_ref','iE','ID','iTwin','variant_SF', 'A_twin','A_residual','A_detwin', 'variant_detwin_TF'};
T5 = cell2table(cell(0,length(variableNames5)));
T5.Properties.VariableNames = variableNames5;

for iE_ref = 3
    iB_ref = iE_ref + 1;
    variant_ref = combined_past_or_present_variant_cell{iB_ref};
    
    for iE = 4:7
        iB = iE + 1;
        variant_residual = combined_current_variant_cell{iB};
        
        % for each parent grain
        for ii = 1:length(gList)
            ID_current = gList(ii);
            ind_g = find(ID_current==gID_0);
            
            euler_current = [gPhi1_0(ind_g), gPhi_0(ind_g), gPhi2_0(ind_g)];
            [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler_current, [0 0 0], [0 0 0], stressTensor, 'Mg', 'twin');
            twin_SFs = abs_schmid_factor(19:24,2);
            
            inds = ismember(ID_overlap, ID_current);
            
            % find child twin grains (from past_or_present_twin_grain_ID)
            children = unique(past_or_present_twin_grain_ID_cell{iB}(inds));
            children(children==0) = [];
            children(isnan(children)) = [];
            
            % find active variants at iE_ref (=3)
            active_variants = unique(variant_ref(inds));
            active_variants(active_variants==0) = [];
            
            % If a past_or_present_twin_grain is completely detwinned, determine its (major) variant and record it      
            detwin_variants = [];
            for jj = 1:length(children)
                ID_children = children(jj);
                ind_children = ismember(past_or_present_twin_grain_ID_cell{iB}, ID_children);
                
                % ==> do we want to ignore very small grains?
                vol_child = sum(ind_children(:));
                
                labels = past_or_present_twin_grain_label_cell{iB}(ind_children);
                
                g1 = 1; 
                g2 = 2;  
                g3 = 3;   
                g4 = 4; 
                g5 = 5;
                
                % If this child/twin grain contains completely detwin past_or_present twin grain, determine its major variant   
                if ismember(g5, labels) && (vol_child>0)
                    variant_ids = combined_past_or_present_variant_cell{iB}(ind_children);
                    major_variant = mode(variant_ids);
                    
                    detwin_variants = unique([detwin_variants, major_variant]);
                end
            end
            
            % check if each variant is significantly detwinned.  
            for iTwin = 1:6
                if ismember(iTwin, active_variants)
                    inds_variant_ref = ismember(ID_overlap, ID_current) & ismember(variant_ref, iTwin);
                    size_variant_ref = sum(inds_variant_ref(:));
                    
                    inds_variant_residual = ismember(ID_overlap, ID_current) & ismember(variant_residual, iTwin);
                    size_variant_residual = sum(inds_variant_residual(:));
                    
                    labels = past_or_present_twin_grain_label_cell{iB}(inds_variant_ref); % g5 = completely detwin grain
                    g5 = 5;
                    size_detwin_grain = sum(labels==g5);
                    
                    if ismember(iTwin, detwin_variants)
                        variant_detwin_TF = true;
                    else
                        variant_detwin_TF = false;
                    end
                    T5 = [T5; {iE_ref, iE, ID_current, iTwin, twin_SFs(iTwin), size_variant_ref, size_variant_residual, size_detwin_grain, variant_detwin_TF}];

                end
            end
            
        end
    end
end


% [A1] Counts of significantly detwinned (and not) vs. twinSF 
for iE = 4:7
    iB = iE + 1;
    
    edges = -0.5:0.05:0.5;
    ind1 = (T5.iE==iE) & (T5.variant_detwin_TF==0);    % 1: NOT significantly detwinned
    ind2 = (T5.iE==iE) & (T5.variant_detwin_TF==1);
    N1 = histcounts(T5.variant_SF(ind1), edges);
    N2 = histcounts(T5.variant_SF(ind2), edges);
    
    figure; hold on; disableDefaultInteractivity(gca);
    hbar = bar(edges(1:end-1)+0.025, [N1(:), N2(:)], 1, 'stacked');
    
    hbar(1).FaceColor = [0 0 1];
    hbar(2).FaceColor = [1 0 0];
    
    set(gca,'fontsize',12, 'XTick',-0.5:0.1:0.5);
    ylim = get(gca,'ylim');
    set(gca, 'ylim', [0 1.5 * ylim(2)]);
    xlabel('Twin Variant Schmid Factor');
    ylabel('Counts');
    
    yyaxis right;
    set(gca, 'ycolor', 'k','fontsize',16);
    ylabel('Percent (%)');
    plot(-0.475:0.05:0.475, N2./(N1+N2) * 100,'-ko','linewidth',1.5);
    ylim = get(gca,'ylim');
    set(gca, 'ylim', [0 1.5 * ylim(2)]);
    
    title(['load step = ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    legend({'Variants not significantly detwinned', 'Variants significantly detwinned','Percent of significantly detwinned'},'Location','northwest');
    print(fullfile(output_dir,['method-2 detwin variant vs SF iE=',num2str(iE),'.tiff']), '-dtiff');
    
    
    %% plot to check. label the grains with 'significantly detwinned variants'
    map = zeros(size(ID_0));
    inds = past_or_present_twin_grain_label_cell{iB}>0;
    map(inds) = past_or_present_twin_grain_label_cell{iB}(inds);
    inds = current_twin_grain_label_cell{iB}>0;
    map(inds) = current_twin_grain_label_cell{iB}(inds);
    
    [f,a,c] = myplot(x,y, map,  boundary_p_iB_to_1_cell{iB});
    
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
    title(['grain summary, load step = ',num2str(iE), ', \epsilon = ',num2str(strain_ebsd(iB),'%.3f')],'fontweight','norma');
    set(gcf,'position',[80,172,1000,600]);
    
    % label those grains
    ids = T5.ID(ind2);
    label_map_with_ID(x,y,ID_overlap, gcf,ids,'r',12,4);
    print(fullfile(output_dir,['method-2 detwin variant grain iE',num2str(iE),'.tiff']), '-dtiff');
    
    close all;
end


