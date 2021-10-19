% Use Mg4Al_U2, iE=3, grain #75, to study if there is a better voting
% method to determine twin variant?
%
% This seems to be a very good example, because
% (1) The child grain has two twin variants, need to seperate
% (2) The orientation has spread, so in the twin area of one variant, the
% pixel-level identification has different results of two variants.
%
% If one 'child grain' contains >1 variants, we should do further analysis 
% 
% This code is used to study the possibility of using the 'connected
% segment length' map to do a further clean up.

%% copy setup
clear; clc; close all;
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD';
save_dir = [working_dir, '\analysis'];
mkdir(save_dir);
% cd(working_dir);

sample_name = 'Mg4Al_U2';

%% copy part-5: find out variants
po_tolerance_angle = 10; % if child grain has misorientation < po_tolerance_angle with undeformed parent grain, it is considered as having parent orientation
twin_tolerance_angle = 10;  % if child grain has misorientation < twin_tolerance_angle to a potential twin variant, it is identified as that twin variant

iE = 3
d = load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_0.mat']));
gID_0 = d.gID;
gPhi1_0 = d.gPhi1;
gPhi_0 = d.gPhi;
gPhi2_0 = d.gPhi2;
ID_0 = d.ID;
boundary_0 = find_boundary_from_ID_matrix(ID_0);

d = load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
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
boundary_c = find_boundary_from_ID_matrix(ID_c);

ID_variant_grain_wise = zeros(size(ID_p));
ID_variant_point_wise = zeros(size(ID_p));
Misorientation_point_wise = zeros(size(ID_p));
ID_variant_by_csl = zeros(size(ID_p));

for ii = 1:length(gID_p)       % 68
    id_p = gID_p(ii);
    disp(['current ID = ',num2str(id_p)]);
    
    
    
    ind_local = ismember(ID_p, id_p); %ismember(ID, [ID_current,ID_neighbor]);
    
    % Make it one data point wider on each side
    indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
    indC_max = min(size(ID_p,2), find(sum(ind_local, 1), 1, 'last')+1);
    indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
    indR_max = min(size(ID_p,1), find(sum(ind_local, 2), 1, 'last')+1);
    
    ID_c_local = ID_c(indR_min:indR_max, indC_min:indC_max);
    x_local = x(indR_min:indR_max, indC_min:indC_max);
    y_local = y(indR_min:indR_max, indC_min:indC_max);
    boundary_c_local = boundary_c(indR_min:indR_max, indC_min:indC_max);
    phi1_c_local = phi1_c(indR_min:indR_max, indC_min:indC_max);
    phi_c_local = phi_c(indR_min:indR_max, indC_min:indC_max);
    phi2_c_local = phi2_c(indR_min:indR_max, indC_min:indC_max);
    
    vg_local = zeros(size(ID_c_local)); % variant map, grain level, local  
    vp_local = zeros(size(ID_c_local)); % variant map, pixel level, local
    
    
        
    % Find the parent orientation at iE = 0
    ind_0 = find(gID_0 == id_p);
    euler_0 = [gPhi1_0(ind_0), gPhi_0(ind_0), gPhi2_0(ind_0)]
    
    if ~isempty(ind_0)
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
        tbl = table(col_ID, col_euler, misorientation); % for debugging
        
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
        for jj = 1:length(id_c)     % 3
            id = id_c(jj);  % twin grains id
            ind = (gID_c == id);
            euler_id = [gPhi1_c(ind), gPhi_c(ind), gPhi2_c(ind)];
            
            
            variant_child_grain = zeros(size(ID_c_local)); % child grain variant map, pixel level, local
            miso_child_grain = zeros(size(ID_c_local)); 
            
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
                
                ind_list_local = find(ID_c_local==id);
                ind_list = find(ID_c==id);
                for mm = 1:length(ind_list)
                    ind = ind_list(mm);
                    euler_c = [phi1_c(ind), phi_c(ind), phi2_c(ind)];
                    
                    misorientation = [];
                    % compare [euler of this pixel]  vs  [euler of the iTwin system of the parent orientation]
                    for kk = 1:6
                        euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                        misorientation(kk) = calculate_misorientation_euler_d(euler_c, euler_kk, 'HCP');
                    end
                    [miso, iVariant] = min(abs(misorientation));
                    Misorientation_point_wise(ind) = miso;
                    ID_variant_point_wise(ind) = iVariant;
                    
                    variant_child_grain(ind_list_local(mm)) = iVariant; 
                    miso_child_grain(ind_list_local(mm)) = miso;
                end
            else
                warning('Twin grain misorientation with parent orientation > 10 deg, rejected as a variant:');
                str = sprintf('iE = %d, ii = %d, ID_parent = %d, jj = %d, ID_twin = %d \n', iE, ii, id_p, jj, id);
                disp(str);
            end            
            active_variants = ismember(1:6, unique(variant_child_grain(:)));
            
            % now we have map_to_analyze = variant_child_grain>0, 'active_variants', 'trace_dir'  
            [abs_schmid_factor, ~, ~] = trace_analysis_TiMgAl(euler_po, [0,0,0], [0,0,0], [-1 0 0; 0 0 0; 0 0 0], 'Mg', 'twin');
            traceDir = abs_schmid_factor(19:24,3); 
            csl_v_map = calculate_csl(variant_child_grain>0, active_variants, traceDir);
            
            vp_local = vp_local + csl_v_map; % add 'child grain variant map by csl' to 'variant pixel local'    
        end
    end
    map = ID_variant_by_csl(indR_min:indR_max, indC_min:indC_max);  % crop map
    map = map + vp_local;   % add grain's variant map
    ID_variant_by_csl(indR_min:indR_max, indC_min:indC_max) = map;  % past back
end

myplot(ID_variant_by_csl,boundary_p);

save(fullfile(save_dir,'try_csl.mat'), 'ID_variant_by_csl');
%% 
myplot(variant_child_grain);
caxis([0 6]);
myplot(miso_child_grain);
myplot(csl_v_map);

% variant_grain_wise{iE} = ID_variant_grain_wise;
% variant_point_wise{iE} = ID_variant_point_wise;




