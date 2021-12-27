% Study the twin contributed strain (corrected by twin SF) vs.
% (o) global strain estimated from EBSD
% (o) twin area fraction

% Method for calculate twin contributed strain:
% twin_contributed_strain = 0.129 * sum_over_grains *
% sum_over_variants(num_of_variant_pixels * SF) / total_num_of_twin_pixels

% Note that the processed data follows the OIM convention for grain
% orientation.  For a grain mostly twinned, the grain orientation is the
% twin orientation.  Therefore, we always need to load data of iE=0 as a
% reference for the parent orientation!!!

addChenFunction;
output_dir = 'E:\zhec umich Drive\0_temp_output\2021-12-22 analyze grains and twins';
mkdir(output_dir);

% [data_dir, sample_name, sample_ID, plot_symbol, group_number]
cells = {'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1;
    'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1;
    'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2;
    'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2;
    'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3;
    'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3;
    'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'o', 4;
    'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 's', 4;
    'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 'x', 4;
    'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'o', 5;
    'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 's', 5;
    'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 'x', 5};

%%
for ii = 1:size(cells,1)
    
    sample_dir = cells{ii,1};
    sample_name = cells{ii,2};
    sample_ID = cells{ii,3};
    
    % load variant maps at all load steps
    d = matfile(fullfile(sample_dir, 'variant_maps.mat'));
    variant_map_cell = d.variant_point_wise;
    % twin data
    d = matfile(fullfile(sample_dir, 'twin_pct.mat'));
    tbl_1 = d.tbl;
    strain_ebsd = tbl_1.strain_ebsd;
    
    % table to summarize twin area, corrected twin area, at each iE
    variable_names = {'iE', 'total_twin_pixels', 'total_twin_pixels_corrected', 'total_pixels', 'twin_pct_corrected', 'twin_strain_corrected'};
    tbl_3 = cell2table(cell(0,length(variable_names)));
    tbl_3.Properties.VariableNames = variable_names;
    
    % reference parent orientation data from iE=0
    d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_0.mat']));
    ID = d.ID;
    gPhi1 = d.gPhi1;
    gPhi = d.gPhi;
    gPhi2 = d.gPhi2;
    gID = d.gID;
    gDiameter = d.gDiameter;
    gArea = d.gArea;
    
    uniqueID = [];    
    % First need to go through all iEs, find unique parent grains to summarize
    for iE = 0:13
        % parent grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        gID_p = d.gID;
        if iE == 0 
            uniqueID = gID_p;
        else
            uniqueID = intersect(uniqueID, gID_p);
        end
    end
    
    % analyze each load step iE
    for iE = 0:13
        iB = iE + 1;
        
        variant_map = variant_map_cell{iB};
        
        % create table, summarize twin by twin
        variable_names = {'iE','ID','iTwin','tSF', 'n_variant_pixels', 'variant_area_corrected', 'gDiameter','gArea'};
        tbl_2 = cell2table(cell(0,length(variable_names)));
        tbl_2.Properties.VariableNames = variable_names;
        
        % parent grain data
        d = matfile(fullfile(sample_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
        ID_p = d.ID;
        % boundary_p = find_one_boundary_from_ID_matrix(ID_p);
        
        % child grain data
        % d = matfile(fullfile(sample_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));
        
        for jj = 1:length(uniqueID)
            ID_current = uniqueID(jj);
            
            ind = find(gID==ID_current);
            
            euler_po = [gPhi1(ind), gPhi(ind), gPhi2(ind)];
            gd = gDiameter(ind);
            ga = gArea(ind);
            
            % note that the euler angles were already corrected.
            [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler_po, [0 0 0], [0 0 0], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
            tSFs = abs_schmid_factor(19:24,2);
            for iTwin = 1:6
                tSF = tSFs(iTwin);
                n_variant_pixels = sum(sum(ismember(ID_p,ID_current) & ismember(variant_map, iTwin)));
                variant_area_corrected = n_variant_pixels * tSF;
                tbl_2 = [tbl_2; {iE, ID_current, iTwin, tSF, n_variant_pixels, variant_area_corrected, gd, ga}];
            end
            
        end
        
        total_twin_pixels = sum(tbl_2.n_variant_pixels);
        total_twin_pixels_corrected = sum(tbl_2.variant_area_corrected);
        total_pixels = numel(ID(:));
        twin_pct_corrected = total_twin_pixels_corrected/total_pixels;
        twin_strain_corrected = twin_pct_corrected * 0.129;
        
        tbl_3 = [tbl_3; {iE, total_twin_pixels, total_twin_pixels_corrected, total_pixels, twin_pct_corrected, twin_strain_corrected}];
        
        save(fullfile(output_dir, [sample_name,'_iE=',num2str(iE),'.mat']),'tbl_2');
    end
    
    tbl = innerjoin(tbl_1, tbl_3);
    
    save(fullfile(output_dir, [sample_name,'.mat']),'tbl');
end

%% (1) plot: x = macroscopic strain, y = twin contributed (corrected) strain
close all;
figure('Position',[500 500 720 400]);
hold on;
legend_str = [];
colors = [0 0 0 ; mix_color(4)];

for ii = 1:length(cells)
    sample_name = cells{ii,2};
    sample_ID = cells{ii,3};
    
    load(fullfile(output_dir, [sample_name,'.mat']),'tbl');
    
    x = abs(tbl.strain_ebsd);
    y = tbl.twin_strain_corrected;
    format_str = cells{ii,4};
    gn = cells{ii,5};
    plot(x(1:4),y(1:4), format_str, 'LineWidth',2, 'MarkerSize',8, 'color', colors(gn,:) );
    % plot(x(5:8),y(5:8), 'or');
    % plot(x(9:11),y(9:11), 'og');
    % plot(x(12:14),y(12:14), 'ob');
    switch gn
        case {1,4}
            legend_str{ii} = [strrep(sample_ID,'_','\_'), ' fine'];
        case 2
            legend_str{ii} = [strrep(sample_ID,'_','\_'), ' medium'];
        case {3,5}
            legend_str{ii} = [strrep(sample_ID,'_','\_'), ' coarse'];
    end
end
set(gca,'ylim', get(gca,'ylim').*[0 1]);
xlabel('Macroscopic Strain');
ylabel('Twin Contributed Strain');
legend(legend_str, 'location','eastoutside');
set(gca,'fontsize',14);

print(fullfile(output_dir, 'twin strain vs global strain.tiff'),'-dtiff');
%% (2) plot: x = twin area fraction, y = twin contributed (corrected) strain
close all;
figure('Position',[500 500 720 400]);
hold on;
legend_str = [];
colors = [0 0 0 ; mix_color(4)];

for ii = 1:length(cells)
    sample_name = cells{ii,2};
    sample_ID = cells{ii,3};
    
    load(fullfile(output_dir, [sample_name,'.mat']),'tbl');
    
    x = abs(tbl.("twinPct %")/100);
    y = tbl.twin_strain_corrected;
    format_str = cells{ii,4};
    gn = cells{ii,5};
    plot(x(1:14),y(1:14), format_str, 'LineWidth',2, 'MarkerSize',8, 'color', colors(gn,:) );
    % plot(x(5:8),y(5:8), 'or');
    % plot(x(9:11),y(9:11), 'og');
    % plot(x(12:14),y(12:14), 'ob');
    switch gn
        case {1,4}
            legend_str{ii} = [strrep(sample_ID,'_','\_'), ' fine'];
        case 2
            legend_str{ii} = [strrep(sample_ID,'_','\_'), ' medium'];
        case {3,5}
            legend_str{ii} = [strrep(sample_ID,'_','\_'), ' coarse'];
    end
end
set(gca,'ylim', get(gca,'ylim').*[0 1]);
xlabel('Twin Area Fraction');
ylabel('Twin Contributed Strain');
legend(legend_str, 'location','eastoutside');
set(gca,'fontsize',14);

print(fullfile(output_dir, 'twin strain vs twin fractin.tiff'),'-dtiff');

