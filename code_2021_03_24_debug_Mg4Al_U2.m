% Look at how pixelwise variant determination
% Conclusion: we might be able to reduce the threshold, using child grains
% with 5 degree misorientation with undeformed parent grain to calculate
% deformed parent orientation. But the improvement is still limited.
%% setup
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD';
save_dir = [working_dir, '\analysis'];
mkdir(save_dir);
% cd(working_dir);

sample_name = 'Mg4Al_U2';

%% from part-5
for iE = 3
    d = load(fullfile(save_dir, 'step-4', [sample_name,'_parent_grain_file_iE_0.mat']));
    gID_0 = d.gID;
    gPhi1_0 = d.gPhi1;
    gPhi_0 = d.gPhi;
    gPhi2_0 = d.gPhi2;
    ID_0 = d.ID;
    boundary_0 = find_boundary_from_ID_matrix(ID_0);
    
    d = load(fullfile(save_dir, 'step-4', [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    x = d.x;
    y = d.y;
    gID_p = d.gID;
    gPhi1_p = d.gPhi1;
    gPhi_p = d.gPhi;
    gPhi2_p = d.gPhi2;
    ID_p = d.ID;
    boundary_p = find_boundary_from_ID_matrix(ID_p);
        
    % data where twins (children) are individially labeled with IDs
    d = load(fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));
    gID_c = d.gID;
    gPhi1_c = d.gPhi1;
    gPhi_c = d.gPhi;
    gPhi2_c = d.gPhi2;
    ID_c = d.ID;
    phi1_c = d.phi1;
    phi_c = d.phi;
    phi2_c = d.phi2;
    
    ID_variant_grain_wise = zeros(size(ID_p));
    ID_variant_point_wise = zeros(size(ID_p));
    Misorientation_point_wise = zeros(size(ID_p));
    
    % Look at ID=75 (ii=68), and 70 (ii=64)
    for ii = 68 %1: length(gID_p)        
        id_p = gID_p(ii);
        disp(['current ID = ',num2str(id_p)]);
        
        % ========== plot for illustration
        inds = ID_p==id_p;
        indC_min = find(sum(inds,1),1,'first');
        indC_max = find(sum(inds,1),1,'last');
        indR_min = find(sum(inds,2),1,'first');
        indR_max = find(sum(inds,2),1,'last');
        ID_c_local = ID_c(indR_min:indR_max, indC_min:indC_max);
        x_local = x(indR_min:indR_max, indC_min:indC_max);
        y_local = y(indR_min:indR_max, indC_min:indC_max);
        boundary_c_local = find_one_boundary_from_ID_matrix(ID_c_local);
        myplotm(ID_c_local, 'x',x_local, 'y',y_local, 'tf',boundary_c_local);
        
        % Find the parent orientation at iE = 0    
        ind_0 = find(gID_0 == id_p);
        euler_0 = [gPhi1_0(ind_0), gPhi_0(ind_0), gPhi2_0(ind_0)];
        
        if ~isempty(ind_0)
            % find out all the id numbers on ID_c overlap with the parent grain on ID_p
            id_c = unique(ID_c(ID_p == id_p));
            
            % ======= plot for illustration
            label_map_with_ID(x_local, y_local, ID_c_local, gcf, id_c, 'r', 18, 1);
            
            % Additionally, the child should have > 70% of its pixels overlap with the potential parent grain
            area_c = zeros(size(id_c));
            overlap_c = zeros(size(id_c));
            for kk = 1:length(area_c)
                area_c(kk) = sum(ID_c(:)==id_c(kk));
                overlap_c(kk) = sum((ID_p(:)==id_p)&(ID_c(:)==id_c(kk)));
            end
            overlap_pct_c = overlap_c./area_c;
            id_c(overlap_pct_c<0.75) = [];
            
            % Modify all the child orientation to crystallographically equivalent one that is closest to euler_0  
            misorientation = [];
            col_ID = [];
            col_euler = [];
            for jj = 1:length(id_c)
                id = id_c(jj);
                ind = (gID_c == id);
                euler_id = [gPhi1_c(ind), gPhi_c(ind), gPhi2_c(ind)];
                
                % =================> Important, here modified the grain info of the child grain      
                euler_id = find_closest_orientation_hcp(euler_id, euler_0);
                gPhi1_c(ind) = euler_id(1);
                gPhi_c(ind) = euler_id(2);
                gPhi2_c(ind) = euler_id(3);
                
                misorientation(jj,1) = calculate_misorientation_euler_d(euler_0, euler_id, 'hcp');
                col_ID = [col_ID; id];
                col_euler = [col_euler; euler_id];
            end
            tbl = table(col_ID, col_euler, misorientation) % for debugging  
            
            % (1) Find out which 'children' is actually the 'parent', by checking the misorientation 
            % Note: multiple 'grains' can have the parent orientation
            inds = find(misorientation < 15);   
            if ~isempty(inds)
                % Use the average of all the 'child grains' with parent orientation.
                id_po = id_c(inds);
                inds_po = find(ismember(gID_c,id_po));
                euler_po = calculate_average_dominant_euler_hcp([gPhi1_c(inds_po), gPhi_c(inds_po), gPhi2_c(inds_po)]);
                    
                % remove the 'child grains' with the 'parent orientation', leave only the true 'twin children'
                id_c(inds) = [];
            else
                % the child maybe fully twinned. Just use euler_0 as euler_po  
                euler_po = euler_0;
                warning('The grain might have completely twinned');
                str = sprintf('iE = %d, ii = %d, ID_parent = %d \n', iE, ii, id_p);
                disp(str);
            end            
                
            % (2) Find out if the remaining child grain is a twin 
            % look at jj=3, child id=356
            jj = find(id_c==356);
            for jj = jj %1:length(id_c)
                id = id_c(jj);  % twin grains id
                ind = (gID_c == id);
                euler_id = [gPhi1_c(ind), gPhi_c(ind), gPhi2_c(ind)];
                
                misorientation = [];
                %  if (demoTF==1)&&(jj==1)&&(ii==ii_demo)
                %       hcp_cell('euler',euler_id,'material','Mg','plotPlane',0,'plotBurgers',0,'plotTrace',0);
                %  end
                
                % Determine if this child grain is a twin variant:  
                % Compare [euler of this child] vs [euler of the iTwin system of the parent orientation]  
                for kk = 1:6
                    euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                    misorientation(kk) = calculate_misorientation_euler_d(euler_id, euler_kk, 'HCP');
                    % if (demoTF==1)&&(jj==1)&&(ii==ii_demo)
                    % 	hcp_cell('euler',euler_po,'material','Mg','ss',24+kk);
                    % 	hcp_cell('euler',euler_kk,'material','Mg','ss',24+kk,'plotPlane',0,'plotBurgers',0,'plotTrace',0);
                    % end
                end
                [min_val, iVariant_child] = min(abs(misorientation));    % find out
                
                % ==============> The child grain may be a twin area containing multiple variants. Assume the child orientation represents at least one true twin orientation  
                % Ff small enough, the child grain should be a twin. Do point-wise analysis   
                if min_val < 10
                    
                    ID_variant_grain_wise(ID_c == id) = iVariant_child; % grain-wise variant map 
                    
                    ind_list = find(ID_c==id);
                    
                    % =========================================   use this to record data for each pixel
                    for k=1:6
                        mis_vi{k} = zeros(size(ID_c_local)); 
                    end
                    % % %
                    
                    for k = 1:length(ind_list)
                        ind = ind_list(k);
                        
                        % convert global ind to local sub 
                        [ir,ic] = ind2sub(size(ID_c),ind);
                        ir = ir - indR_min + 1;
                        ic = ic - indC_min + 1;
                    
                        euler_c = [phi1_c(ind), phi_c(ind), phi2_c(ind)]; 
                        
                        misorientation = [];
                        % compare [euler of this pixel]  vs  [euler of the iTwin system of the parent orientation]  
                        for kk = 1:6
                            euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                            misorientation(kk) = calculate_misorientation_euler_d(euler_c, euler_kk, 'HCP');
                            mis_vi{kk}(ir,ic) = misorientation(kk);
                        end
                        [miso, iVariant] = min(abs(misorientation));
                        Misorientation_point_wise(ind) = miso;
                        ID_variant_point_wise(ind) = iVariant;

                    end                    
                else
                    warning('Twin grain misorientation with parent orientation > 10 deg, rejected as a variant:');
                    str = sprintf('iE = %d, ii = %d, ID_parent = %d, jj = %d, ID_twin = %d \n', iE, ii, id_p, jj, id);
                    disp(str);
                end
            end
            
            myplot(mis_vi{3}); caxis([0 5]); title('misorientation for variant 3');
            myplot(mis_vi{6}); caxis([0 5]); title('misorientation for variant 6');
            myplot(ID_variant_point_wise(indR_min:indR_max, indC_min:indC_max));
        end
    end
    
    
    
%     variant_grain_wise{iE} = ID_variant_grain_wise;
%     variant_point_wise{iE} = ID_variant_point_wise;
%     
%     % After processing this iE, Update the modified grain data and save
%     save_dir_5 = [save_dir,'\step-5'];
%     mkdir(save_dir_5);
%     copyfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), ...
%             fullfile(save_dir_5, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
%     gPhi1 = gPhi1_c;
%     gPhi = gPhi_c;
%     gPhi2 = gPhi2_c;
%     save(fullfile(save_dir_5, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), 'gPhi1','gPhi','gPhi2', '-append');
%     
% 
%     myplot(x,y, ID_variant_grain_wise, boundary_p); caxis([0 6]);
%     title(['iE=',num2str(iE)]);
%     print(fullfile(save_dir,['variant_grain_wise_iE=',num2str(iE),'.tif']),'-dtiff');
%     close;
%     
%     myplot(x,y, ID_variant_point_wise, boundary_p); caxis([0 6]);
%     title(['iE=',num2str(iE)]);
%     print(fullfile(save_dir,['variant_pt_wise_iE=',num2str(iE),'.tif']),'-dtiff');
%     close;      

end

% save(fullfile(save_dir,'variant_maps.mat'),'variant_grain_wise','variant_point_wise');
% 
% % copy the updated file back to the main folder
% for iE = 1:13
%     copyfile(fullfile(save_dir_5, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), ...
%         fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
% end
