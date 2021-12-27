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

clear; clc;
addChenFunction;
input_dir = 'E:\zhec umich Drive\0_temp_output\2021-12-22 analyze grains and twins';
output_dir = 'E:\zhec umich Drive\0_temp_output\2021-12-22 twin contributed strain';
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

%% (1) plot: x = macroscopic strain, y = twin contributed (corrected) strain
close all;
figure('Position',[500 500 720 400]);
hold on;
legend_str = [];
colors = [0 0 0 ; mix_color(4)];

for ii = 1:length(cells)
    sample_name = cells{ii,2};
    sample_ID = cells{ii,3};
    
    load(fullfile(input_dir, [sample_name,'.mat']),'tbl');
    
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
    
    load(fullfile(input_dir, [sample_name,'.mat']),'tbl');
    
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

