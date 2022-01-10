%% [Fig 6] Pct variant twinned at selected iEs, histcounts at SF bins as background
% Note that the processed data follows the OIM convention for grain
% orientation.  For a grain mostly twinned, the grain orientation is the
% twin orientation.  Therefore, we always need to load data of iE=0 as a
% reference for the parent orientation!!!

% purpose: look at how we should determine if a variant is considered as
% 'twinned'

addChenFunction;
output_dir = 'E:\zhec umich Drive\0_temp_output\2022-01-04 all SF distribution';
mkdir(output_dir);

% [data_dir, sample_name, sample_ID, plot_symbol, group_number, ymax]
cells = {'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1, 300;
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1, 300;
    'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2, 350;
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2, 350;
    'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3, 200;
    'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3, 200;
    % 'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'x', 4, 450;
    'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 'o', 4, 450;
    'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 's', 4, 450;
    % 'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'x', 5, 200;
    'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 'o', 5, 200;
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 's', 5, 200};

variant_map_dir = 'E:\zhec umich Drive\All twin variant maps cleaned';
plot_process = 1;
mkdir(fullfile(output_dir, 'process'));

%% Summarize, but need to consider these problem related to data clean up
% (1) Twin areas that are too small might need to be ignored (those <16
% data points should be due to clean up) --> done
% (2) We need a reasonable method to determine: which variants are active
% in a grain?  If we just look at twinned pixels, they might be isolated
% pixels rather than 'grains'.  They might be noise pixels.
% So, we currently use this solution:
% for each child grain, if >10% of its pixels are variant_i, then active
% (3) polishing twins needs to be corrected.
% However, for variants that exist at load step 0, we need to check these
% twins again. If at load steps 1,2,3, pixes of this variant significantly
% increase (>5% of grain size), these can be considered as real twins, not
% polishing.
% Also, create a 'twin area size map' to show how many pixels there are for
% the twins ?
% (4) We tend to use parent grain at undeformed load step. However, when
% look at grain size, edge grains could increase size due to more area
% visible in EBSD map.  So, maybe use deformed parent grain size?

edges = -0.5:0.05:0.5;
for ii = 1:size(cells,1)
    
    sample_dir = cells{ii,1};
    sample_name = cells{ii,2};
    sample_ID = cells{ii,3};
    
    % load variant maps at all load steps
    % d = matfile(fullfile(sample_dir, 'variant_maps.mat'));
    d = matfile(fullfile(variant_map_dir, [sample_name,'_variant_maps.mat']));
    vp_map_cell = d.variant_point_wise; % variant pixel-level
    vg_map_cell = d.variant_grain_wise; % variant grain-level
    
    % reference parent orientation data from iE=0
    d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_0.mat']));
    ID = d.ID;
    gPhi1 = d.gPhi1;
    gPhi = d.gPhi;
    gPhi2 = d.gPhi2;
    gID = d.gID;
    gDiameter = d.gDiameter;
    gArea = d.gArea;
    
    pg_list = [];   % parent grain list
    % First need to go through all iEs, find unique parent grains to summarize
    for iE = 0:13
        % parent grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        gID_p = d.gID;
        if iE == 0
            pg_list = gID_p;
        else
            pg_list = intersect(pg_list, gID_p);
        end
    end
    
    pretwin_id_itwin_list = [0, 0]; % initialize
    % analyze each load step iE
    for iE = 0:13
        disp([sample_name,' iE=',num2str(iE)]);
        iB = iE + 1;
        
        vp_map = vp_map_cell{iB};
        vg_map = vg_map_cell{iB};
        
        variant_pct_map = zeros(size(vp_map));
        
        [nR,nC] = size(vp_map);
        
        % create table, summarize twin by twin
        variable_names = {'iE','ID','iTwin','tSF', 'n_variant_pixels', 'variant_area_corrected', 'gDiameter','gArea'};
        tbl_2 = cell2table(cell(0,length(variable_names)));
        tbl_2.Properties.VariableNames = variable_names;
        
        % parent grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        boundary_p = find_one_boundary_from_ID_matrix(ID_p);
        gID_p = d.gID;
        gDiameter_p = d.gDiameter;
        gArea_p = d.gArea;
        % determine resolution
        y = d.y;
        um_per_dp = y(2) - y(1);
    
        % child grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));
        ID_c = d.ID;
        
        % for each parent grain
        for ip = 1:length(pg_list)
            ID_current = pg_list(ip);
            
            ind = find(gID==ID_current);
            euler_po = [gPhi1(ind), gPhi(ind), gPhi2(ind)];
            
            % for grain size and grain diameter, use deformed parent grain might be better    
            ind_p = find(gID_p==ID_current);
            gd = gDiameter_p(ind_p);
            ga = gArea_p(ind_p);
            
            % find child grains
            ind_p = ismember(ID_p,ID_current);
            cg_list = unique(ID_c(ind_p));
            
            iv_activeTF = zeros(1,6); % variant active TF
            for ic = 1:length(cg_list)
                id_c = cg_list(ic);
                ind_c = ismember(ID_c, id_c);
                
                cg_size = sum(ind_c(:)); % child grain size
                for iTwin = 1:6
                    ind_v = ismember(ID_c, id_c)&ismember(vp_map, iTwin);
                    iv_size = sum(ind_v(:));    % i variant size
                    variant_pct_map(ind_v) = iv_size/cg_size;
                    % If larger than 10% of child grain is this variant
                    if iv_size/cg_size > 0.1
                        iv_activeTF(iTwin) = 1;
                    end
                end
            end
            % remove the pre-existing ones
            inds = pretwin_id_itwin_list(:,1) == ID_current;
            pre_twins = pretwin_id_itwin_list(inds,2);
            iv_activeTF(pre_twins) = 0;
            
            
            % note that the euler angles were already corrected.
            [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler_po, [0 0 0], [0 0 0], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
            tSFs = abs_schmid_factor(19:24,2);
            for iTwin = 1:6
                tSF = tSFs(iTwin);
                if iv_activeTF(iTwin)==1
                    inds = ismember(ID_p,ID_current) & ismember(vp_map, iTwin);
                    n_variant_pixels = sum(inds(:));
                else
                    n_variant_pixels = 0;
                end
                variant_area_corrected = n_variant_pixels * tSF;
                tbl_2 = [tbl_2; {iE, ID_current, iTwin, tSF, n_variant_pixels, variant_area_corrected, gd, ga}];
            end
        end
        
        % if iE==0, check and record pre-existing twins, then remove all from table
        if iE == 0
            if plot_process==1
                % print what it is like before remove
                ind = (tbl_2.iE==iE)&(tbl_2.n_variant_pixels>0);    % twinned variants
                ind2 = (tbl_2.iE==iE)&(tbl_2.n_variant_pixels==0);
                [N_t,~] = histcounts(tbl_2.tSF(ind), edges);
                [N_nt,~] = histcounts(tbl_2.tSF(ind2), edges);
                
                figure; hold on; disableDefaultInteractivity(gca);
                hbar = bar(edges(1:end-1)+0.025, [N_nt(:), N_t(:)], 1, 'stacked');
                
                ymax = cells{ii,6};
                set(gca,'fontsize',12, 'XTick',-0.5:0.1:0.5,'ylim',[0 ymax]);
                
                xlabel('Twin Variant Schmid Factor');
                ylabel('Counts');
                
                yyaxis right;
                set(gca, 'ycolor', 'k','fontsize',16, 'ylim',[0 100]);
                ylabel('Percent (%)');
                plot(-0.475:0.05:0.475, N_t./(N_t+N_nt) * 100,'-ko','linewidth',1.5);
                
                title([sample_ID,', iE = ',num2str(iE)],'fontweight','normal');
                legend({'Variants Not Twinned', 'Variants Twinned','Percent of Variants Twinned'},'Location','northwest');
                
                print(fullfile(output_dir, 'process', [sample_ID,'  iE = ',num2str(iE),'_b.tiff']),'-dtiff');
            end
            
            % Check these twins again. If at load steps 1,2,3, pixes of this
            % variant significantly increase (>5% of grain size), these can
            % be considered as real twins, not polishing-twins.
            tbl_t = tbl_2;  % make a copy, to label real-twins
            tbl_t.n_variant_pixels = zeros(size(tbl_t,1),1);
            for kk = 1:size(tbl_2,1)
                id = tbl_2.ID(kk);
                itwin = tbl_2.iTwin(kk);
                npxls = tbl_2.n_variant_pixels(kk);
                for ie = 1:3
                    ib = ie + 1;
                    % parent grain data (deformed)
                    d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(ie),'.mat']));
                    ID_local = d.ID;
                    vp_local = vp_map_cell{ib};
                    
                    npxls_g = ismember(ID_local, id);   % # of pixels in grain  
                    npxls_g = sum(npxls_g(:));
                    
                    npxls_v = ismember(ID_local, id) & ismember(vp_local, itwin);
                    npxls_v = sum(npxls_v(:));  % # of pixels in twin variant  
                    
                    if (npxls_v - npxls)/npxls_g > 0.05
                        tbl_t.n_variant_pixels(kk) = tbl_2.n_variant_pixels(kk);     % Make the variant volume 0 in the summary, which means pre-existing
                    end
                end
            end
            % find the polishing twins: has volume in tbl_2, but not in tbl_t  
            ind = (tbl_2.n_variant_pixels>0)&(tbl_t.n_variant_pixels==0);
            pretwin_id_itwin_list = [tbl_2.ID(ind), tbl_2.iTwin(ind)];
            
            tbl_2 = tbl_t;  % copy from the modified, leave only real twins              
            
            % plot pre-existing twin variant map
            if plot_process == 1
                ie = 0;
                ib = ie + 1;
                % parent grain data (deformed)
                d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(ie),'.mat']));
                ID_local = d.ID;
                boundary_local = find_one_boundary_from_ID_matrix(ID_local); 
                vp_local = vp_map_cell{ib};
                mask_ok = zeros(size(vp_local));
                for mm = 1:size(pretwin_id_itwin_list,1)
                   id = pretwin_id_itwin_list(mm,1);
                   itwin = pretwin_id_itwin_list(mm,2);
                   ind = ismember(ID_local, id) & ismember(vp_local, itwin);
                   mask_ok(ind) = 1;
                end
                vp_local(~mask_ok) = 0;
                myplot(vp_local, boundary_local); caxis([0 6]);
                print(fullfile(output_dir, 'process', [sample_name, '_pre-existing.tiff']),'-dtiff');
            end
        end
        
        if plot_process == 1
            map_t = vp_map;
            map_t(variant_pct_map>0.1) = 0;
            myplot(map_t,boundary_p); caxis([0 6]);
            print(fullfile(output_dir, 'process', [sample_name,' removed_small_v_iE=',num2str(iE),'.tif']),'-dtiff');
            
            edges = -0.5:0.05:0.5;
            ind = (tbl_2.iE==iE)&(tbl_2.n_variant_pixels>0);    % twinned variants
            ind2 = (tbl_2.iE==iE)&(tbl_2.n_variant_pixels==0);
            [N_t,~] = histcounts(tbl_2.tSF(ind), edges);
            [N_nt,~] = histcounts(tbl_2.tSF(ind2), edges);
            
            figure; hold on; disableDefaultInteractivity(gca);
            hbar = bar(edges(1:end-1)+0.025, [N_nt(:), N_t(:)], 1, 'stacked');
            ymax = cells{ii,6};
            set(gca,'fontsize',12, 'XTick',-0.5:0.1:0.5,'ylim',[0 ymax]);
            xlabel('Twin Variant Schmid Factor');
            ylabel('Counts');
            
            yyaxis right;
            set(gca, 'ycolor', 'k','fontsize',16, 'ylim',[0 100]);
            ylabel('Percent (%)');
            plot(-0.475:0.05:0.475, N_t./(N_t+N_nt) * 100,'-ko','linewidth',1.5);
            
            title([sample_ID,', iE = ',num2str(iE)],'fontweight','normal');
            legend({'Variants Not Twinned', 'Variants Twinned','Percent of Variants Twinned'},'Location','northwest');
            print(fullfile(output_dir, 'process', [sample_ID,'  iE = ',num2str(iE),'.tiff']),'-dtiff');
            
            close all;
        end
        save(fullfile(output_dir, [sample_name,'_iE=',num2str(iE),'.mat']), 'tbl_2', 'pretwin_id_itwin_list');
    end
    
end
close all;

%% [1 plot active variant pct vs SF bins] Load each sample. Then hold on, plot curve for iEs = 0:3, 4:7, 8:10, 11:13
mkdir(fullfile(output_dir, 'active pct'));
close all;
edges = -0.5:0.05:0.5;
cs = inferno(5);
colors = [cs(1:4,:);
    cs(1:4,:);
    cs(1:3,:);
    cs(1:3,:)];

markers = {'-.o','-.d','-.s','-.^', ...
    '-.o','-.d','-.s','-.^', ...
    '-.o','-.d','-.s', ...
    '-.o','-.d','-.s', ...
    '-o','-d','-s','-^'};

for ii=1:length(edges)-1
    xstr{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
end

for ii = 1:size(cells,1)
    
    sample_dir = cells{ii,1};
    sample_name = cells{ii,2};
    sample_ID = cells{ii,3};
    
    clear pct;
    for iE = 0:13
        iB = iE + 1;
        load(fullfile(output_dir, [sample_name,'_iE=',num2str(iE),'.mat']),'tbl_2');
        
        ind = (tbl_2.n_variant_pixels>0);
        ind2 = (tbl_2.n_variant_pixels==0);
        [N_t,~] = histcounts(tbl_2.tSF(ind), edges);
        [N_nt,~] = histcounts(tbl_2.tSF(ind2), edges);
        d_int = (edges(3)-edges(2))/2;
        xpos = [edges(2)-d_int, edges(2:end-1)+d_int];
        
        pct(iB,:) = N_t./(N_t+N_nt);
    end
    
    % just the histogram, using bar plot, show SF distribution
    figure;
    hbar = bar(xpos, [N_nt(:)+N_t(:)], 1);
    set(gca,'ycolor', 'k');
    set(gca, 'fontsize',14, 'XTick', -0.5:0.1:0.5, 'xlim',[-0.3 0.5]);
    hbar(1).FaceColor = [.85 .85 .85];
    xlabel('Twin Variant Schmid Factor');
    ylabel('Counts');
    print(fullfile(output_dir, 'active pct', [sample_name, ' SF distribution.tiff']),'-dtiff');
    
    iEs = {0:3, 4:7, 8:10, 11:13};
    
    for ii = 1:4
        figure; hold on;
        yyaxis right;
        hbar = bar(xpos, [N_nt(:)+N_t(:)], 1);
        set(gca,'ycolor', [.4 .4 .4]);
        set(gca, 'fontsize',13, 'XTick', -0.5:0.1:0.5, 'xlim',[-0.3 0.5]);
        hbar(1).FaceColor = [.85 .85 .85];
        xlabel('Twin Variant Schmid Factor');
        ylabel('Counts');
        
        yyaxis left;
        set(gca, 'ycolor', 'k', 'ylim',[0 100]);
        ylabel('Percent (%)');
        legend_str = [];
        for iE = iEs{ii}
            iB = iE + 1;
            plot(edges(1:end-1)+0.025, 100*pct(iB,:), markers{iB}, 'color', colors(iB,:), 'linewidth', 1.5);
            legend_str = [legend_str, {['Load step ',num2str(iE)]}];
        end
        set(gca, 'SortMethod', 'depth');
        legend(legend_str, 'position',[0.52, 0.72, 0.27, 0.18]);
        print(fullfile(output_dir, 'active pct', [sample_name, ' SF distribution ',num2str(ii),'.tiff']),'-dtiff');
    end
    close all;
end

%% [double check the strange thing] We would like to highlight the negative SF twins in the coarse grain Mg
mkdir(fullfile(output_dir, 'negative SF twin map'));

edges = -0.5:0.05:0.5;
for ii = 1:size(cells,1)
    
    sample_dir = cells{ii,1};
    sample_name = cells{ii,2};
    sample_ID = cells{ii,3};
    
    % load variant maps at all load steps
    % d = matfile(fullfile(sample_dir, 'variant_maps.mat'));
    d = matfile(fullfile(variant_map_dir, [sample_name,'_variant_maps.mat']));
    vp_map_cell = d.variant_point_wise; % variant pixel-level
    vg_map_cell = d.variant_grain_wise; % variant grain-level
    
    % reference parent orientation data from iE=0
    d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_0.mat']));
    ID = d.ID;
    gPhi1 = d.gPhi1;
    gPhi = d.gPhi;
    gPhi2 = d.gPhi2;
    gID = d.gID;
    gDiameter = d.gDiameter;
    gArea = d.gArea;
    
    pg_list = [];   % parent grain list
    % First need to go through all iEs, find unique parent grains to summarize
    for iE = 0:13
        % parent grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        gID_p = d.gID;
        if iE == 0
            pg_list = gID_p;
        else
            pg_list = intersect(pg_list, gID_p);
        end
    end
    
    pretwin_id_itwin_list = [0, 0]; % initialize
    % analyze each load step iE
    for iE = 0:13
        disp([sample_name,' iE=',num2str(iE)]);
        iB = iE + 1;
        
        vp_map = vp_map_cell{iB};
        vg_map = vg_map_cell{iB};
        
        % with variable_names = {'iE','ID','iTwin','tSF', 'n_variant_pixels', 'variant_area_corrected', 'gDiameter','gArea'};
        load(fullfile(output_dir, [sample_name,'_iE=',num2str(iE),'.mat']),'tbl_2', 'pretwin_id_itwin_list');
        
        [nR,nC] = size(vp_map);
        
        % parent grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        boundary_p = find_one_boundary_from_ID_matrix(ID_p);
        
        % child grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));
        ID_c = d.ID;
        
        % for each parent grain
        for ip = 1:length(pg_list)
            ID_current = pg_list(ip);
            
            % pre-existing twins
            pre_twins_TF = zeros(1,6); % variant active TF
            inds = pretwin_id_itwin_list(:,1) == ID_current;
            pre_twins = pretwin_id_itwin_list(inds,2);
            pre_twins_TF(pre_twins) = 1;
            
            for iTwin = 1:6
                indt = (tbl_2.ID==ID_current)&(tbl_2.iTwin==iTwin);
                tSF = tbl_2.tSF(indt);
                
                % if positive twin SF, remove from the 'variant map'
                if (tSF > 0) || pre_twins_TF(iTwin)==1
                    ind_twin = ismember(ID_p, ID_current) & ismember(vp_map, iTwin);
                    vp_map(ind_twin) = 0;
                end
            end
        end
        
        % after checking all grains, plot and save
        myplot(vp_map, boundary_p);
        print(fullfile(output_dir, 'negative SF twin map', ['negative SF twin ', sample_ID,'  iE = ',num2str(iE),'.tiff']),'-dtiff');
        close all;
    end
    
end
close all;


%% [2 box plot: twin fraction vs SF] ==> actually, this AREA-based analysis, it might be better to do based on the 'variant_pixel_mode' after twin_evolution_analysis.
mkdir(fullfile(output_dir, 'boxplot'));
close all;
edges = -0.5:0.05:0.5;

for ii = 1:size(cells,1)
    
    sample_dir = cells{ii,1};
    sample_name = cells{ii,2};
    sample_ID = cells{ii,3};
    
    % determine resolution
    d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_0.mat']));
    y = d.y;
    um_per_dp = y(2) - y(1);
    
    for iE = 0:13
        iB = iE + 1;
        
        % load table
        load(fullfile(output_dir, [sample_name,'_iE=',num2str(iE),'.mat']),'tbl_2');
        
        tbl_2.vA = tbl_2.n_variant_pixels * um_per_dp^2;
        if any(tbl_2.vA > tbl_2.gArea)
            disp(find(tbl_2.vA > tbl_2.gArea));
            error(' ');
        end
        
        ind = (tbl_2.n_variant_pixels>0);  % active twin variants
        t2 = tbl_2(ind,:);
        vx = t2.tSF;
        vy = t2.n_variant_pixels * um_per_dp^2 ./ t2.gArea;
        
        ylabel_str = ['Twin Variant Area Fraction', char(10), 'Normalized by Grain Area'];
                
        gv = discretize(vx, edges);
        nGroups = length(edges)-1;
        clear labels;
        for ii = 1:length(edges)-1
            labels{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
        end
        figure;disableDefaultInteractivity(gca);
        boxplot([vy; nan*ones(nGroups, 1)], [gv; (1:nGroups)'],'notch','off');
        xlabel('Twin Variant Schmid Factor'); ylabel(ylabel_str)
        set(gca,'xticklabels',labels,'xticklabelrotation',45,'fontsize',15, 'xlim',[8.5 20.5], 'ylim', [-0.1 1.1]);
        title_str = [sample_ID, ' iE=',num2str(iE)];
        title(title_str,'fontweight','normal');
        print(fullfile(output_dir, 'boxplot', [title_str,'.tiff']), '-dtiff');
        
        close all;
    end
end

%% Plot twin pixel Schmid Factor (especially for sample UM129 C3) 
mkdir(fullfile(output_dir, 'SF map'));
for ii = 1:size(cells,1)
    
    sample_dir = cells{ii,1};
    sample_name = cells{ii,2};
    sample_ID = cells{ii,3};
        
    % variant data
    d = matfile(fullfile(variant_map_dir, [sample_name,'_variant_maps.mat']));
    vp_map_cell = d.variant_point_wise; % variant pixel-level
    vg_map_cell = d.variant_grain_wise;
        
    % determine resolution
    d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_0.mat']));
    y = d.y;
    um_per_dp = y(2) - y(1);
    
    for iE = 0:13
        iB = iE + 1;
        
        vp_map = vp_map_cell{iB};
        
        [nR,nC] = size(y);
        sfMap = nan*zeros(nR,nC);   % initialize
        
        % parent
        d = matfile(fullfile(sample_dir,[sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        boundary_p = find_one_boundary_from_ID_matrix(ID_p);
    
        load(fullfile(output_dir, [sample_name,'_iE=',num2str(iE),'.mat']),'tbl_2');

        ind = (tbl_2.n_variant_pixels>0);  % active twin variants
        t2 = tbl_2(ind,:);
        
        for jj = 1:size(t2,1)
            tSF = t2.tSF(jj);
            tFraction = t2.n_variant_pixels(jj) * um_per_dp^2 ./ t2.gArea(jj);
            iTwin = t2.iTwin(jj);
            ID_current = t2.ID(jj);
            inds = ismember(ID_p, ID_current) & ismember(vp_map, iTwin);
            
            sfMap(inds) = tSF;
        end
        
        myplot(sfMap, boundary_p);  
        title(' ');
        grid off; 
        cmap = parula(256);
        cmap = cmap(1:240,:);
        colormap(cmap);
        caxis([-0.2, 0.5]);
        print(fullfile(output_dir, 'SF map', [sample_name, ' SF map iE_',num2str(iE),'.tiff']), '-dtiff');
        close all;
    end
end



