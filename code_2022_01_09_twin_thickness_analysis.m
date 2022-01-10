% analyze twin thickness

clear; clc; close all;
addChenFunction;

% [data_dir, sample_name, sample_ID, plot_symbol, group_number, ymax, sample_material]
cells = {'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1, 50, 'Mg4Al FG';
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1, 50, 'Mg4Al FG';
    'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2, 350, 'Mg4Al MG';
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2, 350, 'Mg4Al MG';
    'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3, 100, 'Mg4Al CG';
    'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3, 100, 'Mg4Al CG';
    % 'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'x', 4, 30, 'Mg FG';
    'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 'o', 4, 30, 'Mg FG';
    'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 's', 4, 30, 'Mg FG';
    % 'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'x', 5, 200, 'Mg CG';
    'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 'o', 5, 200, 'Mg CG';
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 's', 5, 200, 'Mg CG'};

% use to determine if a child grain is a twin grain
variant_map_dir = 'E:\zhec umich Drive\All twin variant maps cleaned';
% use to load pretwin_id_itwin_list for polishing induced twins
polishing_twin_dir = 'E:\zhec umich Drive\0_temp_output\2022-01-04 all SF distribution';

% location of the twin evolution (twin detwin retwin) data
input_dir = 'E:\zhec umich Drive\0_temp_output\all twin evolution analysis';

output_dir = 'E:\zhec umich Drive\0_temp_output\twin thickness analysis';
mkdir(output_dir);
mkdir(fullfile(output_dir, 'maps'));

iE_max = 3; % we don't always need all load steps
%%
close all;

for icell = 1:size(cells,1)
    
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    sample_ID = cells{icell,3};
    
    % load variant maps at all load steps
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
    
    variable_names = {'iE','ID','id_c','iTwin','gx','gy','a','b','alpha','tSF','gd','ga'};
    tbl = cell2table(cell(0,length(variable_names)));
    tbl.Properties.VariableNames = variable_names;
    
    
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
    
    
    for iE = 0:iE_max
        iB = iE + 1;
        
        % load polishing induced [id,itwin] list
        d = matfile(fullfile(polishing_twin_dir, [sample_name,'_iE=',num2str(iE),'.mat']));
        pretwin_id_itwin_list = d.pretwin_id_itwin_list;
        
        % we need parent grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        boundary_p = find_one_boundary_from_ID_matrix(ID_p);
        gID_p = d.gID;
        gDiameter_p = d.gDiameter;
        gArea_p = d.gArea;
        y_c = d.y;    % determine resolution
        um_per_dp = y_c(2) - y_c(1);
        
        % these are the child grain data
        data = grain_file_to_data(fullfile(sample_dir, '..', [sample_name,' grain_file_type_1 iE=',num2str(iE),'.txt']), ...
            fullfile(sample_dir, '..', [sample_name,' grain_file_type_2 iE=',num2str(iE),'.txt']));
        
        ID_c = data.ID;
        x_c = data.x;
        y_c = data.y;
        
        gID_c = data.gID;
        gCenterX_c = data.gCenterX;
        gCenterY_c = data.gCenterY;
        gAspectRatio_c = data.gAspectRatio;
        gMajorAxis_c = data.gMajorAxis;
        gMinorAxis_c = data.gMinorAxis;
        gMajorOrientation_c = data.gMajorOrientation;
        gMaxFeretDia_c = data.gMaxFeretDia;
        gMinFeretDia_c = data.gMinFeretDia;
        
        vp_map = vp_map_cell{iB};
        boundary_c = find_one_boundary_from_ID_matrix(ID_c);
        
        
        % for each parent grain, find all child grain
        % for each child grain, find major iTwin, and make sure it is not polishing twin
        % summarize
        for ip = 1:length(pg_list)
            id_p = pg_list(ip);
            
            ind = find(gID==id_p);
            euler_po = [gPhi1(ind), gPhi(ind), gPhi2(ind)];
            [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler_po, [0 0 0], [0 0 0], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
            tSFs = abs_schmid_factor(19:24,2);
            
            % for grain size and grain diameter, use deformed parent grain might be better
            ind_p = find(gID_p==id_p);   % parent grain id, deformed
            gd = gDiameter_p(ind_p);
            ga = gArea_p(ind_p);
            
            % find child grains
            ind_p = ismember(ID_p, id_p);
            cg_list = unique(ID_c(ind_p));
            
            % find variant active TF
            iv_activeTF = zeros(1,6);
            for ic = 1:length(cg_list)
                id_c = cg_list(ic);
                ind_c = ismember(ID_c, id_c);
                
                iTwin = mode(vp_map(ind_c)); % mode of child grain's twin variant label
                
                % if iTwin not in pretwin list, then OK
                if iTwin>0 && ~ismember([id_p,iTwin], pretwin_id_itwin_list, 'rows')
                    ind = find(gID_c == id_c);
                    gx = gCenterX_c(ind);
                    gy = gCenterY_c(ind);
                    
                    a = gMajorAxis_c(ind);
                    b = gMinorAxis_c(ind);
                    alpha = gMajorOrientation_c(ind);
                    tsf = tSFs(iTwin);
                    
                    % append to table
                    tbl = [tbl; {iE, id_p, id_c, iTwin, gx, gy, a, b, alpha, tsf, gd, ga}];
                end
            end
        end
        
        % plot to illustrate
        myplot(x_c,y_c,vp_map,boundary_c); caxis([0 6]);
        make_variant_map_background('bg',[1,1,1],'tickOff',false);
        xlabel('x (\mum)');
        ylabel('y (\mum)');
        title(['load step ',num2str(iE)], 'fontweight','normal');
        set(gca,'fontsize',14);
        
        ind = tbl.iE==iE;
        tbl_local = tbl(ind,:);
        for ir = 1:size(tbl_local,1)
            a = tbl_local.a(ir);
            b = tbl_local.b(ir);
            gx = tbl_local.gx(ir);
            gy = tbl_local.gy(ir);
            alpha = tbl_local.alpha(ir);
            plot3_ellipse(a,b,gx,gy,alpha);
        end
        print(fullfile(output_dir, 'maps', [sample_name, '_fitted_twin_iE_',num2str(iE),'.tiff']),'-dtiff');
        close all;
    end
    
    save(fullfile(output_dir, [sample_name, '_twin_thickness_table.mat']), 'tbl');
end

%% [explore 1] Boxplot: absolute_twin_thickness is different in different materials. Compare at given iE
mkdir(fullfile(output_dir,'abs thick vs mater'));
close all;
% plot for each load step
for iE = 1:3
    iB = iE + 1;
    
    thickness = [];
    gN = [];
    gLabel = [];
    
    iG = 1; % record current group number
    
    for icell = [2,5,7,10] %1:size(cells,1)
        sample_dir = cells{icell,1};
        sample_name = cells{icell,2};
        sample_ID = cells{icell,3};
        sample_material = cells{icell,7};
        % load tbl, for sample, for all iEs
        load(fullfile(output_dir, [sample_name, '_twin_thickness_table.mat']), 'tbl');
        
        ind = (tbl.iE==iE);
        thickness = [thickness; tbl.b(ind) * 2];
        gN = [gN; repmat(iG, sum(ind),1)];
        gLabel{iG} = sample_material;
        iG = iG + 1;
    end
    
    figure;
    boxplot(thickness, gN);
    title(['Load Step ',num2str(iE)], 'fontweight', 'normal');
    ylabel('Twin Thickness (\mum)');
    set(gca,'XTickLabel',gLabel, 'XTickLabelRotation',45);
    set(gca,'fontsize',14, 'ylim',[0 quantile(thickness,0.987)]);
    print(fullfile(output_dir,'abs thick vs mater', ['iE_',num2str(iE),'.tiff']),'-dtiff');
end
close all;

%% [explore 2] Boxplot: 2D grian diameter normalized twin thickness. The fine grain alloys show higher normalized thickness.  Compare at given iE.
mkdir(fullfile(output_dir,'normalized thick vs mater'));
close all;
% plot for each load step
for iE = 1:3
    iB = iE + 1;
    
    thickness = [];
    gd = [];
    thick_norm = [];
    gN = [];
    gLabel = [];
    
    iG = 1; % record current group number
    
    for icell = [2,5,7,10] %1:size(cells,1)
        sample_dir = cells{icell,1};
        sample_name = cells{icell,2};
        sample_ID = cells{icell,3};
        sample_material = cells{icell,7};
        % load tbl, for sample, for all iEs
        load(fullfile(output_dir, [sample_name, '_twin_thickness_table.mat']), 'tbl');
        
        ind = (tbl.iE==iE);
        thickness = [thickness; tbl.b(ind) * 2];
        gd = [gd; tbl.gd(ind)];
        thick_norm = [thick_norm; tbl.b(ind) * 2 ./ tbl.gd(ind)];
        gN = [gN; repmat(iG, sum(ind),1)];
        gLabel{iG} = sample_material;
        iG = iG + 1;
    end
    
    figure;
    boxplot(thick_norm, gN);
    title(['Load Step ',num2str(iE)], 'fontweight', 'normal');
    ylabel(['Twin Thickness Normalized',char(10),'by 2D Grain Diameter']);
    set(gca,'XTickLabel',gLabel, 'XTickLabelRotation',45);
    set(gca,'fontsize',14, 'ylim',[0 quantile(thick_norm,0.987)]);
    print(fullfile(output_dir,'normalized thick vs mater', ['iE_',num2str(iE),'.tiff']),'-dtiff');
end
close all;

%% [explore 3] boxplot: [twin thickness] vs [SF]. Mg4Al has more outliers at high SF
mkdir(fullfile(output_dir,'thick vs SF'));
edges = -0.5:0.05:0.5;
for icell =  [2,5,7,10] %1:size(cells,1)
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    sample_ID = cells{icell,3};
    % load tbl, for sample, for all iEs
    load(fullfile(output_dir, [sample_name, '_twin_thickness_table.mat']), 'tbl');
    
    for iE = 1:3
        iB = iE + 1;
        
        ind = (tbl.iE==iE);
        vx = tbl.tSF(ind);
        vy = tbl.b(ind)*2;
        
        gv = discretize(vx, edges);
        nGroups = length(edges)-1;
        ymax = cells{icell,6};
        
        gv = discretize(vx, edges);
        nGroups = length(edges)-1;
        clear labels;
        for ii = 1:length(edges)-1
            labels{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
        end
        figure;disableDefaultInteractivity(gca);
        boxplot([vy; nan*ones(nGroups, 1)], [gv; (1:nGroups)'],'notch','off');
        xlabel('Twin Variant Schmid Factor');
        ylabel('Twin Thickness (\mum)');
        set(gca,'xticklabels',labels,'xticklabelrotation',45,'fontsize',15, 'xlim',[10.5 20.5],'ylim',[0 quantile(vy,0.99)]);
        title_str = [sample_ID, ' iE=',num2str(iE)];
        title(title_str,'fontweight','normal');
        print(fullfile(output_dir,'thick vs SF', [sample_name,'_iE_',num2str(iE),'.tiff']), '-dtiff');
        close all;
    end
end

%% [explore 4] boxplot: [2D grain diameter normalized twin thickness] vs [SF]. Mg4Al has more outliers at high SF.   pure Mg have almost no effect.
mkdir(fullfile(output_dir,'normalized thick vs SF'));
edges = -0.5:0.05:0.5;
for icell =  [2,5,7,10] %1:size(cells,1)
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    sample_ID = cells{icell,3};
    % load tbl, for sample, for all iEs
    load(fullfile(output_dir, [sample_name, '_twin_thickness_table.mat']), 'tbl');
    
    for iE = 1:3
        iB = iE + 1;
        
        ind = (tbl.iE==iE);
        vx = tbl.tSF(ind);
        vy = tbl.b(ind)*2 ./ tbl.gd(ind);
        
        gv = discretize(vx, edges);
        nGroups = length(edges)-1;
        ymax = cells{icell,6};
        
        gv = discretize(vx, edges);
        nGroups = length(edges)-1;
        clear labels;
        for ii = 1:length(edges)-1
            labels{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
        end
        figure;disableDefaultInteractivity(gca);
        boxplot([vy; nan*ones(nGroups, 1)], [gv; (1:nGroups)'],'notch','off');
        xlabel('Twin Variant Schmid Factor');
        ylabel(['Twin Thickness Normalized',char(10),'by 2D Grain Diameter']);
        set(gca,'xticklabels',labels,'xticklabelrotation',45,'fontsize',15, 'xlim',[10.5 20.5],'ylim',[0 0.6]);
        title_str = [sample_ID, ' iE=',num2str(iE)];
        title(title_str,'fontweight','normal');
        print(fullfile(output_dir,'normalized thick vs SF', [sample_name,'_iE_',num2str(iE),'.tiff']), '-dtiff');
        close all;
    end
end

%% [explore 5] # of twins per grain area

mkdir(fullfile(output_dir,'twin density'));
close all;
% plot for each load step

nPerA = [];
gN = [];
gLabel = [];
iG = 1; % record current group number

for icell = [2,5,7,10] %1:size(cells,1)
    sample_dir = cells{icell,1};
    sample_name = cells{icell,2};
    sample_ID = cells{icell,3};
    sample_material = cells{icell,7};
    
    % load tbl, for sample, for all iEs
    load(fullfile(output_dir, [sample_name, '_twin_thickness_table.mat']), 'tbl');
    for iE = 1:3
        iB = iE + 1;
        
        ind = (tbl.iE==iE);
        tbl_local = tbl(ind,:);
        tbl_local = sortrows(tbl_local, {'ID','iTwin'}, {'ascend','ascend'});
        
        % count for each grain, how many twins
        nR = size(tbl_local,1);
        count_col = zeros(nR,1);   % initialize, record total number of child twins in each parent grain. End row = count, other rows = 0.
        % algorithm: for each row: if next row the same, count++, row ++;
        % else record count, row++ and count=1.
        indr = 1;
        count = 1;
        while indr<nR
            if tbl_local.ID(indr+1)==tbl_local.ID(indr) % && tbl_local.iTwin(indr+1)==tbl_local.iTwin(indr)
                count = count + 1;
            else
                count_col(indr) = count;
                count = 1;
            end
            indr = indr + 1;
        end
        count_col(indr) = count; % record last row
        tbl_local.count = count_col;
        
        ind = tbl_local.count>0;
        tbl_local = tbl_local(ind,:);
        
        nPerA = [nPerA; tbl_local.count ./ tbl_local.ga];
        gN = [gN; repmat(iG, size(tbl_local,1), 1)];
        gLabel{iG} = [sample_material, ' Ls',num2str(iE)];
        iG = iG + 1;
    end
    
end

figure;
boxplot(nPerA, gN);
ylabel(['# Twins per 2D Grain Area (1/\mum^2)']);
set(gca,'XTickLabel',gLabel, 'XTickLabelRotation',45, 'yTick',[1e-5,1e-4,1e-3,1e-2]);
set(gca,'fontsize',12,  'yscale','log');
title(' ', 'fontweight', 'normal');

print(fullfile(output_dir,'twin density', ['twin density.tiff']),'-dtiff');




