
%% setup
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2020-12-05 UM_134 Mg_C1 insitu EBSD';
save_dir = [working_dir, '\analysis'];
cd(working_dir);

%% Find out the twin variants 
for iE = 1:13
    [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_1_parent.txt']));
    [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_2_parent.txt']));
    
    % find column
    column_index_1 = find_variable_column_from_grain_file_header(EBSD_header_1,...
        {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
    column_index_2 = find_variable_column_from_grain_file_header(EBSD_header_2,...
        {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});
    
    gID_p = EBSD_data_2(:,column_index_2(1));
    gPhi1_p = EBSD_data_2(:,column_index_2(2));
    gPhi_p = EBSD_data_2(:,column_index_2(3));
    gPhi2_p = EBSD_data_2(:,column_index_2(4));
    gArea_p = EBSD_data_2(:,column_index_2(9));
    
    % EBSD data, from type-1 grain file. (column, data) pair:
    % (1,phi1) (2,phi) (3,phi2) (4,xMicron) (5,yMicron) (6,IQ) (7,CI) (8,Fit) (9,grain-Id) (10,edgeGrain?)
    % Read EBSD data.  IQ,CI,Fit are not needed for now, but might need in future
    x = EBSD_data_1(:,column_index_1(5));
    y = EBSD_data_1(:,column_index_1(6));
    unique_x = unique(x(:));
    ebsdStepSize = unique_x(2) - unique_x(1);
    mResize = (max(x(:)) - min(x(:)))/ebsdStepSize + 1;
    nResize = (max(y(:)) - min(y(:)))/ebsdStepSize + 1;
    
    x = reshape(EBSD_data_1(:,column_index_1(5)),mResize,nResize)';
    y = reshape(EBSD_data_1(:,column_index_1(6)),mResize,nResize)';
    ID_p = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
    boundary_p = find_boundary_from_ID_matrix(ID_p);
    
    % overwrite data
    load(fullfile(save_dir, ['data_with_ID_overwrite_iE_',num2str(iE),'.mat']), 'ID','gID','gPhi1','gPhi','gPhi2','gNeighbors');
    ID_p = ID;
    gID_p = gID;
    gPhi1_p = gPhi1;
    gPhi_p = gPhi;
    gPhi2_p = gPhi2;
    
    % data where twins (children) are individially labeled with IDs
    [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_1.txt']));
    [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,['UM134_Mg_C1_iE=',num2str(iE),' grain_file_type_2.txt']));
    
    gID_c = EBSD_data_2(:,column_index_2(1));
    gPhi1_c = EBSD_data_2(:,column_index_2(2));
    gPhi_c = EBSD_data_2(:,column_index_2(3));
    gPhi2_c = EBSD_data_2(:,column_index_2(4));
    gArea_c = EBSD_data_2(:,column_index_2(9));
    
    ID_c = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
    
    % Example, biggest grain, iE=3, ii=71, id_p = 70, id_c = [133,135 (parent), 136,150,153,162,163,180 (children)]
    % iE=5, id_p = 69, id_c = [144,148 (parent), 146,160,172,179,193 (children)]
    
    demoTF = 0;
    ii_demo = 71;
    % show this grain where (1) the whole grain has an ID, (2) each twin area has individual ID
    if demoTF==1
        myplot(ID_p, boundary_p);
        set(gca,'xlim',[185,345], 'ylim', [80,240],'fontsize',18,'xticklabel',[],'yticklabel',[]);
        caxis([30,90]);
        title('');
        label_map_with_ID(x,y,ID_p,gcf,[70],'k',24, 2);
        
        myplot(ID_c, boundary_p);
        set(gca,'xlim',[185,345], 'ylim', [80,240],'fontsize',18,'xticklabel',[],'yticklabel',[]);
        caxis([132,168]);
        title('');
        label_map_with_ID(x,y,ID_c,gcf,[133,135,136,150,153,162,163,180],'k',24, 2);
    end
    
    ID_variant = zeros(size(ID_p));
    for ii = 1:length(gID_p)
        id_p = gID_p(ii);
        ind_p = find(gID_p == id_p);
        euler_p = [gPhi1_p(ind_p), gPhi_p(ind_p), gPhi2_p(ind_p)];
        
        % find out all the id numbers on ID_c overlap with the parent grain on ID_p
        id_c = unique(ID_c(ID_p == id_p));
        
        % (1) find out which 'children' is actually the 'parent', by checking the misorientation
        if length(id_c) > 1
            misorientation = [];
            col_ID = [];
            col_euler = [];
            for jj = 1:length(id_c)
                id_twin = id_c(jj);
                ind_twin = (gID_c == id_twin);
                euler_id = [gPhi1_c(ind_twin), gPhi_c(ind_twin), gPhi2_c(ind_twin)];
                misorientation(jj) = calculate_misorientation_euler_d(euler_p, euler_id, 'hcp');
                col_ID = [col_ID; id_twin];
                col_euler = [col_euler; euler_id];
            end
            col_misorientation = misorientation';
            tbl = table(col_ID,col_euler,col_misorientation);
            
            % multiple 'grains' can have the parent orientation
            inds = find(misorientation < 10);
            % use the 1st one for the parent orientation
            id_parent_orientation = id_c(inds(1));  % id on ID_c with the parent orientation
            euler_parent_orientation = [gPhi1_c(id_parent_orientation), gPhi_c(id_parent_orientation), gPhi2_c(id_parent_orientation)];
            
            % [Important modification] use the closest orientation to the parent orientation, without considering symmetry
            euler_parent_orientation = find_closest_orientation_hcp(euler_parent_orientation, euler_p);
            
            % remove the 'twin grains' with the 'parent orientation', leave only the true 'twin children'
            id_c(inds) = [];
            
            % (2) Find out the corresponding twin variant for each 'twin grian'
            for jj = 1:length(id_c)
                id_twin = id_c(jj);  % twin grains id
                ind_twin = (gID_c == id_twin);
                euler_twin = [gPhi1_c(ind_twin), gPhi_c(ind_twin), gPhi2_c(ind_twin)];
                misorientation = [];
                
                if (demoTF==1)&&(jj==1)&&(ii==ii_demo)
                    hcp_cell('euler',euler_twin,'setting',2,'material','Mg','plotPlane',0,'plotBurgers',0,'plotTrace',0);
                end
                % compare euler of the twin grain  vs  euler of the iTwin system of the parent orientation
                for kk = 1:6
                    euler_kk = euler_by_twin(euler_parent_orientation, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                    misorientation(kk) = calculate_misorientation_euler_d(euler_twin, euler_kk, 'HCP');
                    if (demoTF==1)&&(jj==1)&&(ii==ii_demo)
                        hcp_cell('euler',euler_parent_orientation,'setting',2,'material','Mg','ss',24+kk);
                        hcp_cell('euler',euler_kk,'setting',2,'material','Mg','ss',24+kk,'plotPlane',0,'plotBurgers',0,'plotTrace',0);
                    end
                end
                [~, iVariant] = min(abs(misorientation));    % find out
                
                % if small enough
                if misorientation(iVariant) < 10
                    ID_variant(ID_c == id_twin) = iVariant;
                else
                    warning(['warning, ID_parent=',num2str(id_p)]);
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

%% This is corrected from geotrans information
% 13 load steps
strain_corrected = [1, 0.9798, 0.9759, ...
    0.9751, 0.9831, 0.9950, 1.0042, ...
    1.0001, 0.9976, 0.9926, 0.9934, ...
    0.9961, 1.0012, 1.0059] - 1;

% strain from strain gage, iE=0:13 + last unloaded step
strain_sg = [0, -0.0035, -0.0037, ...
    -0.0032, -0.0020, -0.0010, 0.0000, ...
    -0.0010, -0.0020, -0.0033, -0.0035, ...
    -0.0020, -0.0010, -0.0000];
stress = [-2, -67, -76, ...
    -4, 55, 67, 99, ...
    -45, -60, -72, -85, ...
    57, 64, 102];

% use corrected, followed by explicitly copied value, 4 segments in 2 cycles   
strain_corrected = [0, -0.0202, -0.0241, ...
    -0.0249, -0.0169, -0.0050, 0.0042, ...
    0.0001, -0.0024, -0.0074, -0.0066, ...
    -0.0039, 0.0012, 0.0059];   
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
           sub_ID_map = ID([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic);
           
           ii = iE + 1;
           twinPct((ir-1)*nc+ic,ii) = sum(sub_variant_map(:)>0)/sum(sub_ID_map(:)>0);
       end
    end
    
end

%% plot
close all;
tAvg = mean(twinPct);
tStd = std(twinPct);

figure; hold on;
colors = parula(5);

inds = {1:3, 3:7, 7:11, 11:14};
errorbar(strain_corrected(inds{4}), 100*tAvg(inds{4}), 100*tStd(inds{4}), '-', 'color',colors(4,:), 'linewidth',1.5,'markersize',18);
errorbar(strain_corrected(inds{3}), 100*tAvg(inds{3}), 100*tStd(inds{3}), '-', 'color',colors(3,:), 'linewidth',1.5,'markersize',18);
errorbar(strain_corrected(inds{2}), 100*tAvg(inds{2}), 100*tStd(inds{2}), '-', 'color',colors(2,:), 'linewidth',1.5,'markersize',18);
errorbar(strain_corrected(inds{1}), 100*tAvg(inds{1}), 100*tStd(inds{1}), '-', 'color',colors(1,:), 'linewidth',1.5,'markersize',18);


% pos = [strain(:),tAvg(:)*100];
% pos(3,:) = [-0.0065, 2];
% pos(4,:) = [-0.014, 8];
% pos(5,:) = [-0.034, 10];
% pos(6,:) = [-0.0225, 9.7];
% pos(7,:) = [-0.027, 3];
% for iE = 1:13
%     ii = iE + 1;
%     text(pos(ii,1),pos(ii,2), [num2str(100*tAvg(ii),2),'\pm',num2str(100*tStd(ii),2),'%'],'fontsize',16);
% end
set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',[-0.03, 0.01],'ylim',[-2 30],'fontsize',18,'fontweight','normal');
xlabel('Global Strain');
ylabel('Twin Area Percent (%)');



