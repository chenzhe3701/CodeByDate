
%% setup
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD';
save_dir = [working_dir, '\analysis'];
cd(working_dir);

%% Find out the twin variants 
for iE = 1:13
    
    d = load(fullfile(save_dir, ['Mg4Al_C3_parent_grain_file_iE_0.mat']));
    gID_0 = d.gID;
    gPhi1_0 = d.gPhi1;
    gPhi_0 = d.gPhi;
    gPhi2_0 = d.gPhi2;
    ID_0 = d.ID;
    boundary_0 = find_boundary_from_ID_matrix(ID_0);
    
    d = load(fullfile(save_dir, ['Mg4Al_C3_parent_grain_file_iE_',num2str(iE),'.mat']));
    x = d.x;
    y = d.y;
    gID_p = d.gID;
    gPhi1_p = d.gPhi1;
    gPhi_p = d.gPhi;
    gPhi2_p = d.gPhi2;
    ID_p = d.ID;
    boundary_p = find_boundary_from_ID_matrix(ID_p);
        
    % data where twins (children) are individially labeled with IDs
    d = load(fullfile(save_dir, ['Mg4Al_C3_grain_file_iE_',num2str(iE),'.mat']));
    gID_c = d.gID;
    gPhi1_c = d.gPhi1;
    gPhi_c = d.gPhi;
    gPhi2_c = d.gPhi2;
    ID_c = d.ID;

    ID_variant = zeros(size(ID_p));
    for ii = 1: length(gID_p)        
        id_p = gID_p(ii);
        ind_p = find(gID_p == id_p);
        euler_p = [gPhi1_p(ind_p), gPhi_p(ind_p), gPhi2_p(ind_p)];
        
        % Find the parent orientation at iE = 0    ==> Also check the misorientation of the parent grain at the two iEs
        ind_0 = find(gID_0 == id_p);
        euler_0 = [gPhi1_0(ind_0), gPhi_0(ind_0), gPhi2_0(ind_0)];
        
        if ~isempty(ind_0)
            calculate_misorientation_hcp(euler_p, euler_0);

            % find out all the id numbers on ID_c overlap with the parent grain on ID_p
            id_c = unique(ID_c(ID_p == id_p));
            
            % Additionally, the child should have > 70% of its pixels overlap with the potential parent grain
            area_c = zeros(size(id_c));
            overlap_c = zeros(size(id_c));
            for kk = 1:length(area_c)
                area_c(kk) = sum(ID_c(:)==id_c(kk));
                overlap_c(kk) = sum((ID_p(:)==id_p)&(ID_c(:)==id_c(kk)));
            end
            overlap_pct_c = overlap_c./area_c;
            id_c(overlap_pct_c<0.75) = [];
            
            % (1) find out which 'children' is actually the 'parent', by checking the misorientation
            if length(id_c) > 1
                misorientation = [];
                col_ID = [];
                col_euler = [];
                for jj = 1:length(id_c)
                    id_twin = id_c(jj);
                    ind_twin = (gID_c == id_twin);
                    euler_id = [gPhi1_c(ind_twin), gPhi_c(ind_twin), gPhi2_c(ind_twin)];
                    misorientation(jj) = calculate_misorientation_euler_d(euler_0, euler_id, 'hcp');
                    col_ID = [col_ID; id_twin];
                    col_euler = [col_euler; euler_id];
                end
                col_misorientation = misorientation';
                tbl = table(col_ID,col_euler,col_misorientation);
                
                % multiple 'grains' can have the parent orientation
                inds = find(misorientation < 15);
                
                if ~isempty(inds)
                    % [option 1] use the 1st one for the parent orientation
                    id_parent_orientation = id_c(inds(1));  % id on ID_c with the parent orientation
                    inds_po = find(gID_c==id_parent_orientation);
                    euler_parent_orientation = [gPhi1_c(inds_po), gPhi_c(inds_po), gPhi2_c(inds_po)];
                    
                    % ===========>[option 2] try to use the average of all the 'grains' with parent orientation.  
                    % Right now, don't know which option is better
                    id_parent_orientation = id_c(inds);  
                    inds_po = find(ismember(gID_c,id_parent_orientation));
                    euler_parent_orientation = ...
                        calculate_average_dominant_euler_hcp([gPhi1_c(inds_po), gPhi_c(inds_po), gPhi2_c(inds_po)]);
                    
                    
                    % [Important modification] use the closest orientation to the parent orientation, without considering symmetry
                    euler_parent_orientation = find_closest_orientation_hcp(euler_parent_orientation, euler_0);
                    
                    % remove the 'twin grains' with the 'parent orientation', leave only the true 'twin children'
                    id_c(inds) = [];
                    
                    % (2) Find out the corresponding twin variant for each 'twin grian'
                    for jj = 1:length(id_c)
                        id_twin = id_c(jj);  % twin grains id
                        ind_twin = (gID_c == id_twin);
                        euler_twin = [gPhi1_c(ind_twin), gPhi_c(ind_twin), gPhi2_c(ind_twin)];
                        misorientation = [];
                        
                        %                 if (demoTF==1)&&(jj==1)&&(ii==ii_demo)
                        %                     hcp_cell('euler',euler_twin,'setting',2,'material','Mg','plotPlane',0,'plotBurgers',0,'plotTrace',0);
                        %                 end
                        % compare euler of the twin grain  vs  euler of the iTwin system of the parent orientation
                        for kk = 1:6
                            euler_kk = euler_by_twin(euler_parent_orientation, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                            misorientation(kk) = calculate_misorientation_euler_d(euler_twin, euler_kk, 'HCP');
                            %                     if (demoTF==1)&&(jj==1)&&(ii==ii_demo)
                            %                         hcp_cell('euler',euler_parent_orientation,'setting',2,'material','Mg','ss',24+kk);
                            %                         hcp_cell('euler',euler_kk,'setting',2,'material','Mg','ss',24+kk,'plotPlane',0,'plotBurgers',0,'plotTrace',0);
                            %                     end
                        end
                        [~, iVariant] = min(abs(misorientation));    % find out
                        
                        % if small enough
                        if misorientation(iVariant) < 10
                            ID_variant(ID_c == id_twin) = iVariant;
                        else
                            warning('Twin grain misorientation with parent orientation > 10 deg, rejected as a variant:');
                            str = sprintf('iE = %d, ii = %d, ID_parent = %d, jj = %d, ID_twin = %d \n', iE, ii, id_p, jj, id_twin);
                            disp(str);
                        end
                    end
                else
                    warning('No parent orientation found:'); % --> but if no parent found, this could be the actual twin ? 2021-02-03 when writing new code.
                    str = sprintf('iE = %d, ii = %d, ID_parent = %d\n', iE, ii, id_p);
                    disp(str);
                end
                
            end
            
        end
    end
    
    variantMap{iE} = ID_variant;
    
    myplot(x,y, ID_variant, boundary_p); caxis([0 6]);
    title(['iE=',num2str(iE)]);
    print(fullfile(save_dir,['variantMap_iE=',num2str(iE),'.tif']),'-dtiff');
    close;
end

save(fullfile(save_dir,'variant_maps.mat'),'variantMap');

%% study twin area fraction
load(fullfile(save_dir,'variant_maps.mat'),'variantMap');


%% Try point-wise analysis
% for iE = 1:13
    
iE = 1

    d = load(fullfile(save_dir, ['Mg4Al_C3_parent_grain_file_iE_0.mat']));
    gID_0 = d.gID;
    gPhi1_0 = d.gPhi1;
    gPhi_0 = d.gPhi;
    gPhi2_0 = d.gPhi2;
    ID_0 = d.ID;
    boundary_0 = find_boundary_from_ID_matrix(ID_0);
    
    d = load(fullfile(save_dir, ['Mg4Al_C3_parent_grain_file_iE_',num2str(iE),'.mat']));
    x = d.x;
    y = d.y;
    gID_p = d.gID;
    gPhi1_p = d.gPhi1;
    gPhi_p = d.gPhi;
    gPhi2_p = d.gPhi2;
    ID_p = d.ID;
    boundary_p = find_boundary_from_ID_matrix(ID_p);
        
    % data where twins (children) are individially labeled with IDs
    d = load(fullfile(save_dir, ['Mg4Al_C3_grain_file_iE_',num2str(iE),'.mat']));
    gID_c = d.gID;
    gPhi1_c = d.gPhi1;
    gPhi_c = d.gPhi;
    gPhi2_c = d.gPhi2;
    ID_c = d.ID;
    phi1_c = d.phi1;
    phi_c = d.phi;
    phi2_c = d.phi2;

    ID_variant = zeros(size(ID_p));
    ID_miso = zeros(size(ID_p));
    [nR,nC] = size(ID_variant);
    for iR = 1:3:nR
        disp(['iR = ',num2str(iR)]);
        for iC = 1:3:nC
            id_p = ID_p(iR,iC);
            % Find the parent orientation at iE = 0    ==> Also check the misorientation of the parent grain at the two iEs
            ind_0 = find(gID_0 == id_p);
            euler_0 = [gPhi1_0(ind_0), gPhi_0(ind_0), gPhi2_0(ind_0)];
            
            id_c = ID_c(iR,iC);
            ind_c = find(gID_c==id_c);
            euler_c = [phi1_c(iR,iC), phi_c(iR,iC), phi2_c(iR,iC)];
            
            if ~isempty(ind_0)
                misorientation = [];
                % compare euler of the twin grain  vs  euler of the iTwin system of the parent orientation
                for kk = 1:6
                    euler_kk = euler_by_twin(euler_0, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                    misorientation(kk) = calculate_misorientation_euler_d(euler_c, euler_kk, 'HCP');
                end
                [miso, iVariant] = min(abs(misorientation)); 
                ID_miso(iR,iC) = miso;
                ID_variant(iR,iC) = iVariant;
            end
        end 
    end
        
    variantMap{iE} = ID_variant;
    
    myplot(x,y, ID_variant, boundary_p); caxis([0 6]);
    title(['iE=',num2str(iE)]);
    print(fullfile(save_dir,['variantMap_pt_wise_iE=',num2str(iE),'.tif']),'-dtiff');
    close;
% end
save(fullfile(save_dir,'variant_maps_pt_wise.mat'),'variantMap');

%%
ID_variant(ID_miso > 10) = 0;
myplot(x,y, ID_variant, boundary_p); caxis([0 6]);

%% This is corrected from geotrans information
% 13 load steps, Mg4Al_C3

strain_corrected = [0, -0.0108, -0.0218, -0.0300, ...
    -0.0296, -0.0241, -0.0120, -0.0025, ...
    -0.0100, -0.0192, -0.0312, ...
    -0.0225, -0.0124, -0.0036];
    
% strain from image analysis, iE=0:13 + last unloaded step
strain = [0, -0.0075, -0.015, -0.025, ...
    -0.023, -0.017, -0.0075, 0, ...
    -0.0075, -0.015, -0.025, ...
    -0.017, -0.0073, 0];

twinPct = zeros(14,8);

for iE = 1:13
    variant_map = variantMap{iE};
    if iE==1
        [nR,nC] = size(variant_map);
    end
    nr = 3;
    nc = 4;
    for ir=1:nr
       for ic = 1:nc
           sub_variant_map = variant_map([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic);
           sub_ID_map = ID_0([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic); % just for counting number of pixels
           
           ii = iE + 1;
           twinPct((ir-1)*nc+ic,ii) = sum(sub_variant_map(:)>0)/sum(sub_ID_map(:)>0);
       end
    end
    
end

%% [1] plot twin pct vs strain gage strain
close all;
tAvg = mean(twinPct);
tStd = std(twinPct);

figure; hold on;
colors = parula(5);

inds = {1:3, 3:7, 7:11, 11:14};
errorbar(strain(inds{1}), 100*tAvg(inds{1}), 100*tStd(inds{1}), '.-', 'color',colors(1,:), 'linewidth',1.5,'markersize',24);
errorbar(strain(inds{2}), 100*tAvg(inds{2}), 100*tStd(inds{2}), '.-', 'color',colors(2,:), 'linewidth',1.5,'markersize',24);
errorbar(strain(inds{3}), 100*tAvg(inds{3}), 100*tStd(inds{3}), '.-', 'color',colors(3,:), 'linewidth',1.5,'markersize',24);
errorbar(strain(inds{4}), 100*tAvg(inds{4}), 100*tStd(inds{4}), '.-', 'color',colors(4,:), 'linewidth',1.5,'markersize',24);

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',[-0.03, 0.005],'ylim',[-2 50],'fontsize',18,'fontweight','normal');
xlabel('Strain from strain gage');
ylabel('Twin Area Percent (%)');
print(fullfile(save_dir,'twin_pct_vs_sg.tiff'),'-dtiff');
%% [2] plot twin pct vs EBSD strain

tAvg = mean(twinPct);
tStd = std(twinPct);

figure; hold on;
colors = parula(5);

inds = {1:3, 3:7, 7:11, 11:14};
errorbar(strain_corrected(inds{1}), 100*tAvg(inds{1}), 100*tStd(inds{1}), '.-', 'color',colors(1,:), 'linewidth',1.5,'markersize',24);
errorbar(strain_corrected(inds{2}), 100*tAvg(inds{2}), 100*tStd(inds{2}), '.-', 'color',colors(2,:), 'linewidth',1.5,'markersize',24);
errorbar(strain_corrected(inds{3}), 100*tAvg(inds{3}), 100*tStd(inds{3}), '.-', 'color',colors(3,:), 'linewidth',1.5,'markersize',24);
errorbar(strain_corrected(inds{4}), 100*tAvg(inds{4}), 100*tStd(inds{4}), '.-', 'color',colors(4,:), 'linewidth',1.5,'markersize',24);

% text(pos(ii,1),pos(ii,2), [num2str(100*tAvg(ii),2),'\pm',num2str(100*tStd(ii),2),'%'],'fontsize',16);

set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',[-0.035, 0.005],'ylim',[-2 50],'fontsize',18,'fontweight','normal');
xlabel('Strain from EBSD');
ylabel('Twin Area Percent (%)');
print(fullfile(save_dir,'twin_pct_vs_ebsd_strain.tiff'),'-dtiff');

