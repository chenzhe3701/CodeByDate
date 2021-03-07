% code_2021_01_31 look at the accuracy of affine transformation estiamted
% strain.
% Method: 
% Based on DIC data, deform EBSD map.
% Use the DIC deformed EBSD map to estimate strain, then compare with DIC
% strain and strain gage strain.

save_dir = 'E:\zhec umich Drive\2021-01-31 evaluate affine';

strain_sg = [0, -0.0075, -0.015, -0.025, -0.023, -0.017, 0];
strain_dic = [0, -0.0012, -0.0117, -0.0186, -0.0172, -0.0124, 0.006];
strain_ebsd = [0, -0.0003, -0.0114, -0.0183, -0.0169, -0.0126, 0.0001];     % this is what the result will be after running this code.

%% load dic, load aligned_EBSD, deform it for each iE using the DIC data at that iE
dic_dir = 'E:\Mg4Al_S1_insitu\SEM Data\stitched_DIC';

dxX = [];
for iE = 0:6    
    disp(iE);
    % load EBSD ID map, which is already aligned to DIC
    load('E:\Mg4Al_S1_insitu\Summary\Mg4Al_S1_EBSD_organized.mat','ID')
    % resize to about 600 x 600
    [nr,nc] = size(ID);
    r = min(floor(nr/600), floor(nc/600));  % ratio needed to resize map to approximately 600 x 600
    ID = ID(1:r:end, 1:r:end);
    
    % boundary = find_one_boundary_from_ID_matrix(ID);
    % myplot(ID, boundary);

   load(fullfile(dic_dir, ['_',num2str(iE),'.mat']), 'x','y','u','v');
   
   x = x(1:r:end, 1:r:end);
   y = y(1:r:end, 1:r:end);
   u = u(1:r:end, 1:r:end);
   v = v(1:r:end, 1:r:end);
    
   % ID deformed to these [x_p, y_p] positions
   x_p = x + u;
   y_p = y + v;
   x_p = inpaint_nans(x_p);
   y_p = inpaint_nans(y_p);
   
   % mid of deformed map
   x_mid = round(nanmean(x_p(:)));
   y_mid = round(nanmean(y_p(:)));
   
   % # rows, cols, and step_size of requested map
   [nR,nC] = size(ID);
   xy_step = x(1,2) - x(1,1);
   
   % Analyze DIC: for each row, find the first non-nan u/v to estimate the
   % elongation, then strain
   elongation = [];
   gage = [];
   for iR = 1:nR
      ind_1 = find(~isnan(u(iR,:)), 1, 'first');
      ind_2 = find(~isnan(u(iR,:)), 1, 'last');
      gage = [gage; x(iR,ind_1), x(iR,ind_2)];
      elongation = [elongation; x(iR,ind_1) + u(iR,ind_1), x(iR,ind_2) + u(iR,ind_2)];
   end
   dx = mean(gage(:,2)) - mean(gage(:,1));
   dX = mean(elongation(:,2)) - mean(elongation(:,1));
   strain_dic_edge(iE+1) = dX/dx - 1;
   dxX = [dxX; dx, dX];
   
   % meshgrid to inquire the deformed ID map
   xi = x_mid - round(nC/2) * xy_step;
   xf = xi + xy_step*(nC-1);
   yi = y_mid - round(nR/2) * xy_step;
   yf = yi + xy_step*(nR-1);
   [x_deformed_mesh, y_deformed_mesh] = meshgrid(xi:xy_step:xf, yi:xy_step:yf);
   
   % reconstructed deformed ID map.
   ID_deformed = interp_data(x_p, y_p, ID, x_deformed_mesh, y_deformed_mesh, [], 'fit', 'nearest');
   % boundary_to = find_one_boundary_from_ID_matrix(ID_to);
   % myplot(ID_to,boundary_to);
   
   
   % rename to regular name, and summarize
   ID = ID_deformed;
   [x,y] = meshgrid(1:nC, 1:nR);
   % summarize gCenter and gEdge
   ids_1 = unique(ID); % all ids
   ids_2 = [unique(ID(1,:)), unique(ID(end,:))];
   ids_3 = [unique(ID(:,1)); unique(ID(:,end))];
   ids_4 = unique([ids_2(:);ids_3(:)]);    % ids on the border
   ids_interior = ids_1(~ismember(ids_1,ids_4));    % adjusted ids that are not on the edge
   
   [gID, ~, IC] = unique(ID);  % IC: each of gID_more_gb's index in ID_more_gb
   accumNumPts = accumarray(IC,1);
   accumX = accumarray(IC,x(:));
   accumY = accumarray(IC,y(:));
   
   gCenterX = accumX./accumNumPts; % reculculate grain center directly from map
   gCenterY = accumY./accumNumPts;
 
   gEdge = ones(size(gID));
   ind = ismember(gID, ids_interior);      % index indicating interior grains
   gEdge(ind) = 0;
   
   save(fullfile(save_dir, ['Mg4Al_S1_artificial_EBSD_iE_',num2str(iE),'.mat']), 'ID', 'x', 'y', 'gID', 'gEdge', 'gCenterX', 'gCenterY');
end

save(fullfile(save_dir, 'Mg4Al_S1_dic_edge_estimated_strain.mat'),'dxX','strain_dic_edge');

%% check the deformed maps. Try use four grains to estimate
IDs = [68, 72;
    223, 230];
strain_four_points = 1;
for iE = 0:6
    load(fullfile(save_dir, ['Mg4Al_S1_artificial_EBSD_iE_',num2str(iE),'.mat']));
    
    gb = find_one_boundary_from_ID_matrix(ID);
    
    myplot(x,y,ID*nan,gb);
    title(['iE=',num2str(iE)]);
    hold on;
    ind = gEdge==1;
    plot(gCenterX(ind), gCenterY(ind), 'ob');
    plot(gCenterX(~ind), gCenterY(~ind), 'or');
    print(fullfile(save_dir,['deformed iE=',num2str(iE)]),'-dtiff');
    close;
    if iE==0
        dx = (gCenterX(gID==IDs(1,2)) - gCenterX(gID==IDs(1,1)) + gCenterX(gID==IDs(2,2)) - gCenterX(gID==IDs(2,1)))/2;
    else
        dX = (gCenterX(gID==IDs(1,2)) - gCenterX(gID==IDs(1,1)) + gCenterX(gID==IDs(2,2)) - gCenterX(gID==IDs(2,1)))/2;
        strain_four_points = [strain_four_points; dX/dx];
    end
end
strain_four_points = strain_four_points - 1
strain_four_points = [0, -0.0002, -0.0133, -0.0199, -0.0188, -0.0138, 0];   % the result

%% use all the non-edge grains to perform affine transformatoin, and display the transformation matrix

for iE = 0:6
    disp(['iE=',num2str(iE)]);
    % load the modified EBSD data at iE
    load(fullfile(save_dir, 'Mg4Al_S1_artificial_EBSD_iE_0.mat'));
    
    gCenterX_0 = gCenterX;
    gCenterY_0 = gCenterY;
    gEdge_0 = gEdge;
    gID_0 = gID;
    
	load(fullfile(save_dir, ['Mg4Al_S1_artificial_EBSD_iE_',num2str(iE),'.mat']));
    
    gID_1 = gID_0(gEdge_0==0);  % interior IDs on the reference map
    gID_2 = gID(gEdge==0);      % interior IDs on the deformed map

    gID_for_affine = intersect(gID_1, gID_2);
    
    loc_0 = ismember(gID_0, gID_for_affine);
    loc_iE = ismember(gID, gID_for_affine);
    
    cpFrom = [gCenterX_0(loc_0),gCenterY_0(loc_0)]; % cpFrom is from deformed image (iE>0)
    cpTo = [gCenterX(loc_iE),gCenterY(loc_iE)];   % cpTo is from ref image (iE=0)
    
    % [2nd] fine tform.
    tform = fitgeotrans(cpFrom, cpTo, 'affine');    % [x_ref, y_ref, 1] * tform.T = [x_def, y_def, 1]
    tforms{iE+1} = tform;     % save the tform
    
end
for iE = 0:6
   tforms{iE+1}.T' 
end


%% This is the result run from the python script
strain_ebsd = [1.0000, 0.9997, 0.9886, 0.9817, 0.9831, 0.9874, 1.0001] - 1








