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

%% [1] code piece from part-5. Look at ID=75 (ii=68). 
% In this grain, the twin area from lower-left to upper-right should be
% variant 6, but many pixels were identified to be variant 3. 
% In this analysis, we try to see if we reduce the tolerance for a child
% grain to be considered as having the parent orientation, how the
% identification results could alter.
% Conclusion: it can reduce the error, but cannot completely solve this
% problem.
po_tolerance_angle = 10; % if child grain has misorientation < po_tolerance_angle with undeformed parent grain, it is considered as having parent orientation
twin_tolerance_angle = 10;  % if child grain has misorientation < twin_tolerance_angle to a potential twin variant, it is identified as that twin variant  
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
    
    % Look at ID=75 (ii=68)
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
            inds = find(misorientation < po_tolerance_angle);   
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
                if min_val < twin_tolerance_angle
                    
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
    
end

%% [2] Now we look at iE=3, grain 30 (ii=29)
% Looks like the variants are more likely to be 1 and 4
% But lots of pixels are identified as variants 2 and 5
% We want to look at how orientation and Schmid factor of these variants
% differ

po_tolerance_angle = 10; % if child grain has misorientation < po_tolerance_angle with undeformed parent grain, it is considered as having parent orientation
twin_tolerance_angle = 10;  % if child grain has misorientation < twin_tolerance_angle to a potential twin variant, it is identified as that twin variant  
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
    
    % Look at ID=75 (ii=68)
    for ii = 29 %1: length(gID_p)        
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
        phi1_c_local = phi1_c(indR_min:indR_max, indC_min:indC_max);
        phi_c_local = phi_c(indR_min:indR_max, indC_min:indC_max);
        phi2_c_local = phi2_c(indR_min:indR_max, indC_min:indC_max);
        
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
            inds = find(misorientation < po_tolerance_angle);   
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
                
            % =========== look at how variant orientations differ
            hcp_cell('euler',euler_po,'material','Mg','ss',24+kk);
            for kk = 1:6
                euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                hcp_cell('euler',euler_kk,'material','Mg','ss',24+kk,'plotPlane',1,'plotBurgers',1,'plotTrace',1);
            end
                
            % (2) Find out if the remaining child grain is a twin 
            for jj = 1:length(id_c)
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
                if min_val < twin_tolerance_angle
                    
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
            
            %%
            [xx,yy] = meshgrid(1:size(ID_c_local,2),1:size(ID_c_local,1));
            [f,a,c] = myplotm(ID_variant_point_wise(indR_min:indR_max, indC_min:indC_max), ...
                'x', xx, 'y', yy, 'tf', boundary_c_local); caxism([0 6]);
            set(gca,'fontsize',16);
            set(c,'limits',[0 6]);
            
            % ============ for selected pixel, check.  pt 1
            indr = 20;
            indc = 40;
            plot(indc,indr,'.r','markersize',32);
            text(indc,indr,'  pt1','fontsize',18, 'color', 'r');
            euler_pt = [phi1_c_local(indr,indc), phi_c_local(indr,indc), phi2_c_local(indr,indc)];
            % check the difference between SF, misorientation
            for kk = 1:6
                euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                misorientation(kk) = calculate_misorientation_euler_d(euler_pt, euler_kk, 'HCP');
            end
            [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler_po, [0 0 0], [0 0 0], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
            t = array2table([(1:6)',abs_schmid_factor(19:24,2),abs_schmid_factor(19:24,3),misorientation(:)]);
            t.Properties.VariableNames = {'vNum','SF','traceDir','misorientation'};
            display(t);
            
            % pt 2
            indc = 25;
            indr = 60;
            plot(indc,indr,'.r','markersize',32);
            text(indc,indr,'  pt2','fontsize',18, 'color', 'r');
            euler_pt = [phi1_c_local(indr,indc), phi_c_local(indr,indc), phi2_c_local(indr,indc)];
            % check the difference between SF, misorientation
            for kk = 1:6
                euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                misorientation(kk) = calculate_misorientation_euler_d(euler_pt, euler_kk, 'HCP');
            end
            [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler_po, [0 0 0], [0 0 0], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
            t = array2table([(1:6)',abs_schmid_factor(19:24,2),abs_schmid_factor(19:24,3),misorientation(:)]);
            t.Properties.VariableNames = {'vNum','SF','traceDir','misorientation'};
            display(t);
            
            % pt 3
            indc = 18;
            indr = 13;
            plot(indc,indr,'.r','markersize',32);
            text(indc,indr,'  pt3','fontsize',18, 'color', 'r');
            euler_pt = [phi1_c_local(indr,indc), phi_c_local(indr,indc), phi2_c_local(indr,indc)];
            % check the difference between SF, misorientation
            for kk = 1:6
                euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                misorientation(kk) = calculate_misorientation_euler_d(euler_pt, euler_kk, 'HCP');
            end
            [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler_po, [0 0 0], [0 0 0], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
            t = array2table([(1:6)',abs_schmid_factor(19:24,2),abs_schmid_factor(19:24,3),misorientation(:)]);
            t.Properties.VariableNames = {'vNum','SF','traceDir','misorientation'};
            display(t);
            
            % pt 4
            indc = 18;
            indr = 18;
            plot(indc,indr,'.r','markersize',32);
            text(indc,indr,'  pt4','fontsize',18, 'color', 'r');
            euler_pt = [phi1_c_local(indr,indc), phi_c_local(indr,indc), phi2_c_local(indr,indc)];
            % check the difference between SF, misorientation
            for kk = 1:6
                euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                misorientation(kk) = calculate_misorientation_euler_d(euler_pt, euler_kk, 'HCP');
            end
            [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler_po, [0 0 0], [0 0 0], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
            t = array2table([(1:6)',abs_schmid_factor(19:24,2),abs_schmid_factor(19:24,3),misorientation(:)]);
            t.Properties.VariableNames = {'vNum','SF','traceDir','misorientation'};
            display(t);
            
        end
    end
    
end















