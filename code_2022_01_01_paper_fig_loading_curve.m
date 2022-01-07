%% [Fig 2] all loading curves

% Suggest to run step by step
addChenFunction;
output_dir = 'E:\zhec umich Drive\0_temp_output\2022 gs alloy paper';
mkdir(output_dir);
clear colors;
colors = [0 0 0; mix_color(4)];

p = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu curve';
f = 'Mg4Al_C1_processed_loading_data.mat';
% d = matfile(fullfile(p,f));
% stress = d.stress;
% strain = d.strain;
% plot(strain, stress, '-', 'color', 'k');
% legend_str = [legend_str, 'Mg4Al C1 fine grain'];

%% location of all data for ref
p = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 insitu curve\borrow some data';
f = 'Mg4Al_U2_processed_loading_data';

p = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu curve';
f = 'Mg4Al_C3_processed_loading_data.mat';

p = 'E:\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu curve';
f = 'Mg4Al_A1_processed_loading_data.mat';

p = 'E:\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu curve';
f = 'Mg4Al_A2_processed_loading_data.mat';

p = 'E:\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu curve';
f = 'Mg4Al_B1_processed_loading_data';

p = 'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu curve';
f = 'Mg4Al_B2_processed_loading_data';

% p = 'E:\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu curve';
% f = 'UM134_Mg_C1_processed_loading_data.mat';

p = 'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu curve';
f = 'UM134_Mg_C2_processed_loading_data.mat';

p = 'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu curve';
f = 'UM134_Mg_C3_processed_loading_data';

% p = 'E:\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu curve';
% f = 'UM129_Mg_C1_processed_loading_data.mat';

p = 'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu curve';
f = 'UM129_Mg_C2_processed_loading_data.mat';

p = 'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu curve';
f = 'UM129_Mg_C3_processed_loading_data.mat';

%%
for n = 1:6
    switch n
        case 1
            % Mg4Al, fine grain
            close all;
            p = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu curve';
            f = 'Mg4Al_C3_processed_loading_data.mat';
            d = matfile(fullfile(p,f));
            stress = d.stress;
            strain = d.strain;
            ind_stop = d.ind_stop;
            title_str = 'Mg4Al C3 fine grain';
            ylim = [-150 150];
        case 2
            % Mg4Al, coarse grain
            p = 'E:\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu curve';
            f = 'Mg4Al_B2_processed_loading_data';
            d = matfile(fullfile(p,f));
            stress = d.stress;
            strain = d.strain;
            ind_stop = d.ind_stop;
            title_str = 'Mg4Al B2 coarse grain';
            ylim = [-150 150];
            
        case 3
            % Pure Mg, fine grain
            p = 'E:\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu curve';
            f = 'UM134_Mg_C3_processed_loading_data';
            d = matfile(fullfile(p,f));
            stress = d.stress;
            strain = d.strain;
            ind_stop = d.ind_stop;
            title_str = 'Mg UM134 C3 fine grain';
            ylim = [-75 75];
            
        case 4
            % Pure Mg, coarse grain
            p = 'E:\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu curve';
            f = 'UM129_Mg_C2_processed_loading_data.mat';
            d = matfile(fullfile(p,f));
            stress = d.stress;
            strain = d.strain;
            ind_stop = d.ind_stop;
            title_str = 'Mg UM129 C2 coarse grain';
            ylim = [-75 75];
            
        case 5
            % Pure Mg, fine grain
            p = 'E:\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu curve';
            f = 'UM134_Mg_C2_processed_loading_data.mat';
            d = matfile(fullfile(p,f));
            stress = d.stress;
            strain = d.strain;
            ind_stop = d.ind_stop;
            title_str = 'Mg UM134 C2 fine grain';
            ylim = [-75 75];
            
        case 6
            % Pure Mg, coarse grain
            p = 'E:\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu curve';
            f = 'UM129_Mg_C3_processed_loading_data.mat';
            d = matfile(fullfile(p,f));
            stress = d.stress;
            strain = d.strain;
            ind_stop = d.ind_stop;
            title_str = 'Mg UM129 C3 coarse grain';
            ylim = [-75 75];
    end
    
    % Normalize stress and strain, calculate the range, so that we plot data
    % point every approximately range/100 distance.
    stress_range = range(stress);
    strain_range = range(strain);
    v_norm = [strain_range, stress_range];
    dmax_norm = pdist2([strain(ind_stop(4)), stress(ind_stop(4))]./v_norm, ...
        [strain(ind_stop(8)), stress(ind_stop(8))]./v_norm, 'euclidean');
    
    close all;
    % cycle-1
    figure; hold on;
    c1 = 1;
    c2 = 4;
    plot(nan,nan,'.','color',colors(c1,:),'markersize',18);
    plot(nan,nan,'.','color',colors(c2,:),'markersize',18);
    
    plot(strain(1), stress(1), 'o', 'color',colors(c1,:), 'linewidth',2);
    
    for ii = 1:7
        ind1 = ind_stop(ii);
        ind2 = ind_stop(ii+1);        
        
        plot(strain(ind1), stress(ind1), '.', 'color',colors(c1,:));
        plot(strain(ind2), stress(ind2), 'o', 'color',colors(c1,:), 'linewidth',2);
        
        pos_old = [strain(ind1), stress(ind1)];
        jj = ind1;
        while jj < ind2
            pos = [strain(jj), stress(jj)];
            if pdist2(pos_old./v_norm, pos./v_norm, 'euclidean') > dmax_norm/100
                pos_old = pos;
                plot(strain(jj), stress(jj), '.', 'color',colors(c1,:), 'markersize', 10);
            end
            jj = jj + 1;
        end
    end
    
    % cycle-2
    for ii = 8:13
        ind1 = ind_stop(ii);
        ind2 = ind_stop(ii+1);
        
        plot(strain(ind1), stress(ind1), '.', 'color',colors(c2,:));
        plot(strain(ind2), stress(ind2), 'o', 'color',colors(c2,:), 'linewidth',2);
        
        pos_old = [strain(ind1), stress(ind1)];
        jj = ind1;
        while jj < ind2
            pos = [strain(jj), stress(jj)];
            if pdist2(pos_old./v_norm, pos./v_norm, 'euclidean') > dmax_norm/100
                pos_old = pos;
                plot(strain(jj), stress(jj), '.', 'color',colors(c2,:), 'markersize', 10);
            end
            jj = jj + 1;
        end
    end
    
    set(gca,'fontsize',16);
    set(gca,'xlim',[-0.03 0.005],'ylim',ylim);
    xlabel('Strain'); ylabel('Stress (MPa)');
    legend('Cycle 1', 'Cycle 2','location','northwest');
    title(title_str, 'fontweight','normal'); %title(' ');
    print(fullfile(output_dir, ['fig 2 ', title_str, '.tif']), '-dtiff');
    
end

close all;