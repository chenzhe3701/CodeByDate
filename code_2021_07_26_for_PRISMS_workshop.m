
% prepare Figures for in-situ EBSD paper and PRISMS workshop 2021
clear; clc; close all;
addChenFunction;
output_dir = 'E:\zhec umich Drive\0_temp_output\for 2021 workshop to edit';
mkdir(output_dir);

%% orientation for the hcp_cell for illustration
close all;
hcp_cell('euler',[1.7 109.5 5], 'plotPlane',1,'plotBurgers',0,'plotTrace',0, 'ss', 30);

%% Use Mg4Al_U2 as example
working_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD';
save_dir = [working_dir, '\analysis'];
sample_name = 'Mg4Al_U2';

%% (task 1) Show IPF map before clean up, and after clean up. Just need simple illustratin
close all;
iE = 1;

% directly read and crop map
% (1) uncleaned map
I = imread(fullfile(working_dir, ['Mg4Al_U2 IPF raw ND iE=',num2str(iE),'.tif']));
figure;
imshow(I);

iR_min = 221;
iR_max = iR_min + 199;
iC_min = 241;
iC_max = iC_min + 199;

I = I(iR_min:iR_max, iC_min:iC_max, :);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, 'fig 1a IPF before clean.tiff'));

% (2) cleaned map
I = imread(fullfile(working_dir, ['Mg4Al_U2 IPF ND iE=',num2str(iE),'.tif']));
figure;
imshow(I);

I = I(iR_min:iR_max, iC_min:iC_max, :);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, 'fig 1b IPF cleaned.tiff'));

%% (task 0) show that the parent/dauther map may be wrong
close all;

iE_a = 1;
iE_b = 3;

% [] ref IPF map at iE=0
I = imread(fullfile(working_dir, 'Mg4Al_U2 parent IPF ND iE=0.tif'));
I = I(iR_min+10 : iR_max+10, iC_min:iC_max, :);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, 'fig 0a IPF iE=0.tiff'));

% [] IPF map at iE_a
I = imread(fullfile(working_dir, ['Mg4Al_U2 parent IPF ND iE=',num2str(iE_a),'.tif']));
I = I(iR_min:iR_max, iC_min:iC_max, :);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, ['fig 0b IPF iE=',num2str(iE_a),'.tiff']));

% [] OIM parent/twin map at iE=2
I = imread(fullfile(working_dir, ['Mg4Al_U2 TwinMap iE=',num2str(iE_a),'.tif']));
I = I(iR_min:iR_max, iC_min:iC_max, :);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, ['fig 0c twin map iE=',num2str(iE_a),'.tiff']));

% [] IPF map at iE_b
I = imread(fullfile(working_dir, ['Mg4Al_U2 parent IPF ND iE=',num2str(iE_b),'.tif']));
I = I(iR_min:iR_max, iC_min:iC_max, :);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, ['fig 0d IPF iE=',num2str(iE_b),'.tiff']));

% [] OIM parent/twin map at iE=3
I = imread(fullfile(working_dir, ['Mg4Al_U2 TwinMap iE=',num2str(iE_b),'.tif']));
I = I(iR_min:iR_max, iC_min:iC_max, :);
figure;
imshow(I);
imwrite(I, fullfile(output_dir, ['fig 0e twin map iE=',num2str(iE_b),'.tiff']));

%% (task 2) show 2 types of grain data: [type-A] ignore twin boundary and only show parent, [type-B] show individual twins. 
iE = 1;

iR_min = 261;
iR_max = iR_min + 99;
iC_min = 241;
iC_max = iC_min + 99;

close all;
% (1) ignore twin boundary, only parent 
d = grain_file_to_data(fullfile(working_dir, ['Mg4Al_U2 parent_grain_file_type_1 iE=',num2str(iE),'.txt']), ...
    fullfile(working_dir, ['Mg4Al_U2 parent_grain_file_type_2 iE=',num2str(iE),'.txt']));

x = d.x(iR_min:iR_max, iC_min:iC_max);
y = d.y(iR_min:iR_max, iC_min:iC_max);
ID = d.ID(iR_min:iR_max, iC_min:iC_max);
boundary = find_one_boundary_from_ID_matrix(ID);
phi1 = d.phi1(iR_min:iR_max, iC_min:iC_max);
phi = d.phi(iR_min:iR_max, iC_min:iC_max);
phi2 = d.phi2(iR_min:iR_max, iC_min:iC_max);

IPF = generate_IPF_map_pixel_wise(x,y, phi1,phi,phi2, boundary, [90,180,0], [0 0 1]);
figure;
imshow(IPF);
% label_map_with_ID(x-min(x(:)),y-min(y(:)),ID,gcf, unique(ID), 'r',18,1);

ID_select = 72;
label_map_with_ID(x-min(x(:)),y-min(y(:)),ID,gcf, ID_select, 'b',18,1);
inds = ismember(ID,ID_select);
set(gcf,'position', [100, 100, 800, 600]);
print(fullfile(output_dir, ['fig 2a ID with label ignore twin iE=',num2str(iE),'.tiff']),'-dtiff','-r300');

imwrite(IPF, fullfile(output_dir, ['fig 2a gID ignore twin iE=',num2str(iE),'.tiff']));

% (2) include twin boundary, show all children 
d = grain_file_to_data(fullfile(working_dir, ['Mg4Al_U2 grain_file_type_1 iE=',num2str(iE),'.txt']), ...
    fullfile(working_dir, ['Mg4Al_U2 grain_file_type_2 iE=',num2str(iE),'.txt']));

x = d.x(iR_min:iR_max, iC_min:iC_max);
y = d.y(iR_min:iR_max, iC_min:iC_max);
ID = d.ID(iR_min:iR_max, iC_min:iC_max);
boundary = find_one_boundary_from_ID_matrix(ID);
phi1 = d.phi1(iR_min:iR_max, iC_min:iC_max);
phi = d.phi(iR_min:iR_max, iC_min:iC_max);
phi2 = d.phi2(iR_min:iR_max, iC_min:iC_max);

IPF = generate_IPF_map_pixel_wise(x,y, phi1,phi,phi2, boundary, [90,180,0], [0 0 1]);
figure;
imshow(IPF);

% IDs_select = unique(ID(inds))
% label_map_with_ID(x-min(x(:)),y-min(y(:)),ID,gcf, IDs_select,'k',18,1);
label_map_with_ID(x-min(x(:)),y-min(y(:)),ID,gcf, 182,'b',18,1);
text(75,37,'177','fontsize',18, 'color','b');
text(82,58,'178','fontsize',18, 'color','b');
set(gcf,'position', [100, 100, 800, 600]);
print(fullfile(output_dir, ['fig 2b ID with label include twin iE=',num2str(iE),'.tiff']),'-dtiff','-r300');

imwrite(IPF, fullfile(output_dir, ['fig 2b gID include twin iE=',num2str(iE),'.tiff']));

%% (task 4) example of: [A] merge grains, [B] divide grains more
% we need to load from processed data

close all;
% (1) show reference load step, parent IPF?, there are 2 parent grains
d = load(fullfile(working_dir, 'analysis', 'step-1', ['Mg4Al_U2_parent_grain_file_iE_0.mat']));
iR_min = 266;
iR_max = iR_min + 169;
iC_min = 226;
iC_max = iC_min + 129;

x = d.x(iR_min:iR_max, iC_min:iC_max);
y = d.y(iR_min:iR_max, iC_min:iC_max);
ID = d.ID(iR_min:iR_max, iC_min:iC_max);
boundary = find_one_boundary_from_ID_matrix(ID);
phi1 = d.phi1(iR_min:iR_max, iC_min:iC_max);
phi = d.phi(iR_min:iR_max, iC_min:iC_max);
phi2 = d.phi2(iR_min:iR_max, iC_min:iC_max);

IPF = generate_IPF_map_pixel_wise(x,y, phi1,phi,phi2, boundary, [90,180,0], [0 0 1]);
figure;
imshow(IPF);
imwrite(IPF, fullfile(output_dir, 'fig 4a ref IPF with gb.tiff'));

%% (2) deformed parent IPF, with grain boundary
iE = 5;
d = load(fullfile(working_dir, 'analysis', 'step-1', ['Mg4Al_U2_parent_grain_file_iE_',num2str(iE),'.mat']));
iR_min = 266;
iR_max = iR_min + 169;
iC_min = 226;
iC_max = iC_min + 129;

x = d.x(iR_min:iR_max, iC_min:iC_max);
y = d.y(iR_min:iR_max, iC_min:iC_max);
ID = d.ID(iR_min:iR_max, iC_min:iC_max);
boundary = find_one_boundary_from_ID_matrix(ID);
phi1 = d.phi1(iR_min:iR_max, iC_min:iC_max);
phi = d.phi(iR_min:iR_max, iC_min:iC_max);
phi2 = d.phi2(iR_min:iR_max, iC_min:iC_max);

IPF = generate_IPF_map_pixel_wise(x,y, phi1,phi,phi2, boundary, [90,180,0], [0 0 1]);
figure;
imshow(IPF);

imwrite(IPF, fullfile(output_dir, ['fig 4b IPF with gb at iE_',num2str(iE),'.tiff']));

%% (3) Simulate the draw mask process, plot: max_misorientation map + mask
tolerance = 3;

ID_local = ID;
x_local = x;
y_local = y;
boundary_local = boundary;

phi1_local = phi1;
phi_local = phi;
phi2_local = phi2;

% for each pixel, find max misorientation within its four neighbors
[nR,nC] = size(ID_local);
misorientation_max = zeros(nR,nC);
for iR = 1:nR
    for iC = 1:nC
        neighbor_orientation = [];
        % add upper neighbor, bottom neighbor, left neighbor, right neighbor
        if iR>1
            neighbor_orientation = [neighbor_orientation; phi1_local(iR-1,iC), phi_local(iR-1,iC), phi2_local(iR-1,iC)];
        end
        if iR<nR
            neighbor_orientation = [neighbor_orientation; phi1_local(iR+1,iC), phi_local(iR+1,iC), phi2_local(iR+1,iC)];
        end
        if iC>1
            neighbor_orientation = [neighbor_orientation; phi1_local(iR,iC-1), phi_local(iR,iC-1), phi2_local(iR,iC-1)];
        end
        if iC<nC
            neighbor_orientation = [neighbor_orientation; phi1_local(iR,iC+1), phi_local(iR,iC+1), phi2_local(iR,iC+1)];
        end
        
        euler_current = [phi1_local(iR,iC), phi_local(iR,iC), phi2_local(iR,iC)];
        
        misorientation = [];
        for ii = 1:size(neighbor_orientation,1)
            misorientation(ii) = calculate_misorientation_euler_d(euler_current, neighbor_orientation(ii,:), 'HCP');
        end
        misorientation_max(iR,iC) = max(misorientation);
    end
end

% hf2 = myplotm(boundary_local,'x',x_local,'y', y_local);
% label_map_with_ID(x_local,y_local,ID_local, gcf, ID_current, 'r');
hf3= myplotm(misorientation_max);
title('');
axis off;
set(gca,'fontsize',18);
caxism([tolerance, 90]);

pos =[
   27.1380   85.5000
   30.5941   78.5879
   62.3896   77.8967
   97.2955   79.2791
  104.2076   84.1176
   91.7658   87.2280
   48.9110   87.2280];
   
h = drawpolygon('Position',pos,'color','r');

print(fullfile(output_dir, 'fig 4c max misorientation map and mask.tiff'), '-dtiff', '-r300');

%% (4) original boundary + new boundary
mask = h.createMask;
boundary_local_new = (misorientation_max > tolerance)&(mask==1);

map = double(boundary);
map(boundary_local_new==1) = 2;
myplotm(map);
caxism([0.1 2.1]);
title('');
axis off;
colorbar off;

print(fullfile(output_dir, 'fig 4d old and new gb map.tiff'), '-dtiff', '-r300');

%% (task 3) Show alignment/correlation of maps at load steps 0 and iE 
iE = 2;
iB = iE + 1;

% load reference data at iE=0
d = load(fullfile(working_dir, 'analysis', 'step-1', ['Mg4Al_U2_parent_grain_file_iE_0.mat']));

gID_0 = d.gID;
gPhi1_0 = d.gPhi1;
gPhi_0 = d.gPhi;
gPhi2_0 = d.gPhi2;
gCenterX_0 = d.gCenterX;
gCenterY_0 = d.gCenterY;
gEdge_0 = d.gEdge;

x = d.x;
y = d.y;
ID_0 = d.ID;
boundary_0 = find_one_boundary_from_ID_matrix(ID_0);

% load deformed iE
d = load(fullfile(working_dir, 'analysis', 'step-1', ['Mg4Al_U2_parent_grain_file_iE_',num2str(iE),'.mat']));

gID = d.gID;
gPhi1 = d.gPhi1;
gPhi = d.gPhi;
gPhi2 = d.gPhi2;
gCenterX = d.gCenterX;
gCenterY = d.gCenterY;
gEdge = d.gEdge;

ID = d.ID;
phi1 = d.phi1;
phi = d.phi;
phi2 = d.phi2;

[boundary, boundaryID, neighborID, tripleTF, tripleID, indTriple, triIDs] = find_one_boundary_from_ID_matrix(ID);

grain_pair{3} = [18, 18;
    31, 29;
    108, 96;
    120, 110];

% (step-2) Rough align using the selected control grains. The result is already decent 
g_0 = grain_pair{iB}(:,1);    % ref at iE=0
g_iE = grain_pair{iB}(:,2);    % iE > 0, considered as deformed

[~, loc_0] = ismember(g_0, gID_0);
[~, loc_iE] = ismember(g_iE, gID);
cpFrom = [gCenterX_0(loc_0),gCenterY_0(loc_0)];     % cpFrom is from ref image (iE=0)
cpTo = [gCenterX(loc_iE),gCenterY(loc_iE)];         % cpTo is from deformed image (iE>0)

% The tform is to transform the [ref @ iE=0] to [deformed iE>0]
tform = fitgeotrans(cpFrom, cpTo, 'affine');    % [x_0, y_0, 1] * tform.T = [x_iE, y_iE, 1]
ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');

% (step-3) Fine align again, based on id_link and use all non-edge grains, and show result! 
[ID_new, id_link_additional, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);

% [[remove]] '0's from linked ids
inds = find(id_link(:,1)==0 | id_link(:,2)==0);
id_link(inds,:) = [];
inds = isnan(ID_0_to_iE)|(ID_new==0); % old ID is nan, or matched ID=0
ID_new(inds) = nan;

% [[remove]] grains that are shown in BOTH the col_2 of id_link and col_2 of in_link_additional  
ind = ismember(id_link(:,2), id_link_additional(:,2));
id_link(ind,:) = [];

% [[remove]] edge grains in id_link 
g_0 = id_link(:,1);
g_iE = id_link(:,2);    
[~, loc_0] = ismember(g_0, gID_0);
[~, loc_iE] = ismember(g_iE, gID);
ind = (gEdge_0(loc_0)==1) | (gEdge(loc_iE)==1);
g_0(ind) = [];
g_iE(ind) = [];

% use cps to transform
[~, loc_0] = ismember(g_0, gID_0);
[~, loc_iE] = ismember(g_iE, gID);
cpFrom = [gCenterX_0(loc_0),gCenterY_0(loc_0)]; % cpFrom is from deformed image (iE>0)
cpTo = [gCenterX(loc_iE),gCenterY(loc_iE)];   % cpTo is from ref image (iE=0)

% The tform is to transform the 'deformed' back into the 'ref'
tform = fitgeotrans(cpFrom, cpTo, 'affine');    % [x_ref, y_ref, 1] * tform.T = [x_def, y_def, 1]

ID_0_to_iE = interp_data(x,y,ID_0, x,y,tform, 'interp', 'nearest');
boundary_0_to_iE = find_boundary_from_ID_matrix(ID_0_to_iE);

% (step-4) try to match grains based on their spatial location, and show the matching result. 
% Transparent grains are the ones without a good match, so the parent grain need to be divided into more grains.     
[ID_new, id_link_additional, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);
% remove 0s
ind = find(id_link(:,1)==0 | id_link(:,2)==0);
id_link(ind,:) = [];
inds = isnan(ID_0_to_iE)|(ID_new==0); % old ID is nan, or matched ID=0
ID_new(inds) = nan;
% ---> just for illustration, I don't want to show the newly assigned ID
inds = ismember(ID_new, id_link_additional(:,1));
ID_new(inds) = nan;

    
% transform maps at iE to i0
boundary_iE_to_0 = interp_data(x,y, boundary, x,y, tform.invert, 'interp', 'nearest');


close all;

% (plot 1, raw) boundary iE=0 black, boundary iE red
figure; hold on;

z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');

z = boundary;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'r');

set(gca,'ydir','reverse', 'fontsize',16);
xlabel('x (\mum)');
ylabel('y (\mum)');
axis equal;
axis square;
print(fullfile(output_dir, ['fig 3a raw overlay iE=0 and ',num2str(iE),'.tiff']), '-dtiff', '-r150');

% (plot 2, aligned) boundary iE=0 black, boundary iE red
figure; hold on;

z = boundary_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'k');

z = boundary_iE_to_0;
z(z==0) = nan;
plot3(x(:),y(:),z(:),'.', 'markersize',1, 'color', 'r');

set(gca,'ydir','reverse', 'fontsize',16);
xlabel('x (\mum)');
ylabel('y (\mum)');
axis equal;
axis square;
print(fullfile(output_dir, ['fig 3b correlated overlay iE=0 and ',num2str(iE),'.tiff']), '-dtiff', '-r150');


%% (task 5) Identification, use parent data in folder 'step-4', child data in folder 'step-2'
close all;
%% (1) show reference load step, parent IPF?, ID_target = 74
ID_target = 74;

d = load(fullfile(working_dir, 'analysis', 'step-4', ['Mg4Al_U2_parent_grain_file_iE_0.mat']));
iR_min = 261;
iR_max = iR_min + 119;
iC_min = 231;
iC_max = iC_min + 119;

x = d.x(iR_min:iR_max, iC_min:iC_max);
y = d.y(iR_min:iR_max, iC_min:iC_max);
ID = d.ID(iR_min:iR_max, iC_min:iC_max);
boundary = find_one_boundary_from_ID_matrix(ID);
phi1 = d.phi1(iR_min:iR_max, iC_min:iC_max);
phi = d.phi(iR_min:iR_max, iC_min:iC_max);
phi2 = d.phi2(iR_min:iR_max, iC_min:iC_max);

IPF = generate_IPF_map_pixel_wise(x,y, phi1,phi,phi2, boundary, [90,180,0], [0 0 1]);
figure;
imshow(IPF);
label_map_with_ID(x-min(x(:)), y-min(y(:)), ID, gcf, ID_target, 'k', 18, 1);
set(gcf,'position',[100 100 800 600]);

print(fullfile(output_dir, 'fig 5a_1 ref IPF iE=0.tiff'), '-dtiff');
imwrite(IPF, fullfile(output_dir, 'fig 5a_2 ref IPF iE=0.tiff'));

gPhi1 = d.gPhi1;
gPhi = d.gPhi;
gPhi2 = d.gPhi2;
gID = d.gID;

ind = (gID==ID_target);
euler_0 = [gPhi1(ind), gPhi(ind), gPhi2(ind)]
hcp_cell('euler',euler_0, 'ss', 1, 'plotPlane',0, 'plotBurgers',0, 'plotTrace',0, 'plotAC',1);

print(fullfile(output_dir, 'fig 5b hcp_cell iE_0.tiff'), '-dtiff');
%% (2) show parent grain at iE = 1
iE = 1;

d = load(fullfile(working_dir, 'analysis', 'step-4', ['Mg4Al_U2_parent_grain_file_iE_',num2str(iE),'.mat']));
iR_min = 256;
iR_max = iR_min + 119;
iC_min = 231;
iC_max = iC_min + 119;

x = d.x(iR_min:iR_max, iC_min:iC_max);
y = d.y(iR_min:iR_max, iC_min:iC_max);
ID = d.ID(iR_min:iR_max, iC_min:iC_max);
boundary = find_one_boundary_from_ID_matrix(ID);
phi1 = d.phi1(iR_min:iR_max, iC_min:iC_max);
phi = d.phi(iR_min:iR_max, iC_min:iC_max);
phi2 = d.phi2(iR_min:iR_max, iC_min:iC_max);

IPF = generate_IPF_map_pixel_wise(x,y, phi1,phi,phi2, boundary, [90,180,0], [0 0 1]);
figure;
imshow(IPF);
label_map_with_ID(x-min(x(:)), y-min(y(:)), ID, gcf, ID_target, 'k', 18, 1);
set(gcf,'position',[100 100 800 600]);

print(fullfile(output_dir, ['fig 5c_1 parent IPF iE=',num2str(iE),'.tiff']), '-dtiff');
imwrite(IPF, fullfile(output_dir, ['fig 5c_2 parent IPF iE=',num2str(iE),'.tiff']));

%% (3) show child grain at iE = 1.  grain IDs = [182,177,178], where 177 is twin

d = load(fullfile(working_dir, 'analysis', 'step-2', ['Mg4Al_U2_grain_file_iE_',num2str(iE),'.mat']));
iR_min = 256;
iR_max = iR_min + 119;
iC_min = 231;
iC_max = iC_min + 119;

x = d.x(iR_min:iR_max, iC_min:iC_max);
y = d.y(iR_min:iR_max, iC_min:iC_max);
ID = d.ID(iR_min:iR_max, iC_min:iC_max);
boundary = find_one_boundary_from_ID_matrix(ID);
phi1 = d.phi1(iR_min:iR_max, iC_min:iC_max);
phi = d.phi(iR_min:iR_max, iC_min:iC_max);
phi2 = d.phi2(iR_min:iR_max, iC_min:iC_max);

IPF = generate_IPF_map_pixel_wise(x,y, phi1,phi,phi2, boundary, [90,180,0], [0 0 1]);
figure;
imshow(IPF);
label_map_with_ID(x-min(x(:)),y-min(y(:)),ID,gcf, 182,'b',18,1);
text(84,50,'177','fontsize',18, 'color','b');
text(93,65,'178','fontsize',18, 'color','b');
set(gcf,'position',[100 100 800 600]);

print(fullfile(output_dir, ['fig 5d_1 child IPF iE=',num2str(iE),'.tiff']), '-dtiff');
imwrite(IPF, fullfile(output_dir, ['fig 5d_2 child IPF iE=',num2str(iE),'.tiff']));

% euler for child grains, id=[182,178] are parent orientation, id=[177] is twin to check  
gPhi1 = d.gPhi1;
gPhi = d.gPhi;
gPhi2 = d.gPhi2;
gID = d.gID;

ind = (gID==182);
euler_a = [gPhi1(ind), gPhi(ind), gPhi2(ind)]
hcp_cell('euler',euler_a, 'ss', 1, 'plotPlane',0, 'plotBurgers',0, 'plotTrace',0);
print(fullfile(output_dir, ['fig 5e hcp_cell iE_',num2str(iE),' ID_182.tiff']), '-dtiff');

ind = (gID==178);
euler_b = [gPhi1(ind), gPhi(ind), gPhi2(ind)]
hcp_cell('euler',euler_b, 'ss', 1, 'plotPlane',0, 'plotBurgers',0, 'plotTrace',0);
print(fullfile(output_dir, ['fig 5f hcp_cell iE_',num2str(iE),' ID_178.tiff']), '-dtiff');

ind = (gID==177);
euler_c = [gPhi1(ind), gPhi(ind), gPhi2(ind)]
hcp_cell('euler',euler_c, 'ss', 1, 'plotPlane',0, 'plotBurgers',0, 'plotTrace',0);
print(fullfile(output_dir, ['fig 5g hcp_cell iE_',num2str(iE),' ID_177.tiff']), '-dtiff');

% averaged parent orientation
euler_po = calculate_average_dominant_euler_hcp([euler_a; euler_b])
euler_po = find_closest_orientation_hcp(euler_po, euler_0)
hcp_cell('euler',euler_po, 'ss', 1, 'plotPlane',0, 'plotBurgers',0, 'plotTrace',0, 'plotAC',1);
print(fullfile(output_dir, ['fig 5h hcp_cell iE_',num2str(iE),'_po.tiff']), '-dtiff');

[S, M] = hcp_symmetry();
for iSym = 1:12
    % equivalent orientations - 2 or 8
    if iSym==2 || iSym==8 || 1
        euler_r = euler_by_M(euler_po, M(:,:,iSym))
        hcp_cell('euler',euler_r, 'ss', 1, 'plotPlane',0, 'plotBurgers',0, 'plotTrace',0, 'plotAC',1);
        print(fullfile(output_dir, ['fig 5h hcp_cell iE_',num2str(iE),'_po eq',num2str(iSym),'.tiff']), '-dtiff');
    end
end

%% (4) show the hcp_cell with avg_parent_orientation, and each of the 6 variant orientation
close all;

% hcp cell of parent orientation, and active TS # 2
hcp_cell('euler',euler_po, 'ss', 26, 'plotPlane',0, 'plotBurgers',0, 'plotTrace',0, 'plotAC',0);
print(fullfile(output_dir, ['fig 5i parent_cell iE_',num2str(iE),'.tiff']), '-dtiff');

% Based on parent orientation, calculate orientation of each possible variant,  
% And calculate misorientation with child grain #177  
misorientation = [];
for kk = 1:6
    euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
    misorientation(kk) = calculate_misorientation_euler_d(euler_c, euler_kk, 'HCP');
    
    hcp_cell('euler',euler_kk, 'ss', 24+kk, 'plotPlane',1, 'plotBurgers',1, 'plotTrace',1, 'plotAC',0);
    print(fullfile(output_dir, ['fig 5i twin_cell iE_',num2str(iE),' ts_',num2str(kk),'.tiff']), '-dtiff');
end

%%  [Task, merge for movie] Mg4Al_U2, IPF map, Twin variant map, Twin type map, curve.  
if 0
    % top modify the curve
    open code_2021_03_01_Mg4Al_U2_curve.m;
    open code_2020_12_29_Mg4Al_C3_curve.m;
    
    % IPF map: directly use tif?
    
    % to modify the twin variant map (not spatially transformed):
    open code_2021_07_30_plot_variant_map_with_bg_color.m
    
    % to modify the twin category map:
    open code_2021_05_04_analyze_persistent_twin.m;
end

% sources of variant map, experimentally measured position
vMap_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis\variant map';
vMap_file = ['variant_pixel_map_iE_',num2str(iE),'.tiff'];

% spatially deformed to iE=0 variant map
vMap_spatially_dir = 'E:\zhec umich Drive\0_temp_output\2021-05-04 Analyze persistent twin\Mg4Al_U2';
vMap_spatially_file = ['pixel variant map iE=',num2str(iE),'.tiff'];

IPF_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD';
IPF_file = ['Mg4Al_U2 parent IPF ND iE=',num2str(iE),'.tif'];

cat_dir = 'E:\zhec umich Drive\0_temp_output\2021-05-04 Analyze persistent twin\Mg4Al_U2';
cat_file = ['frd twinned iE=',num2str(iE),'.tiff'];     % fresh, re, de-twins 

curve_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 insitu curve'; % this sample, lost some loading data
curve_dir = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu curve';
curve_file = 'stress vs strain gage.tiff';

%% For each iE, merge the maps
for iE = 0:13
    close all;
    % (1) Loading curve
    image_1 = imread(fullfile(curve_dir, curve_file));
    % (2) IPF
    IPF_file = ['Mg4Al_U2 parent IPF ND iE=',num2str(iE),'.tif'];
    image_2 = imread(fullfile(IPF_dir, IPF_file));
    % (3) twin variant map
    vMap_file = ['variant_pixel_map_iE_',num2str(iE),'.tiff'];
    image_3 = imread(fullfile(vMap_dir, vMap_file));
    % (4) category map
    cat_file = ['frd twinned iE=',num2str(iE),'.tiff']; 
    image_4 = imread(fullfile(cat_dir, cat_file));
    
%     figure; imshow(image_1);
%     figure; imshow(image_2);
%     figure; imshow(image_3);
%     figure; imshow(image_4);

    I_1 = imcrop(image_1, [0,0, inf,inf]);
    I_1 = imresize(I_1, 620/656);
    
    I_2 = imcrop(image_2, [0,0, inf,inf]);
    I_2 = imresize(I_2, [761, 761]);
    
    I_3 = imcrop(image_3, [200,70, 880,765]);
    I_3 = imresize(I_3, 761/765);
    
    I_4 = imcrop(image_4, [215,75, inf,760]);
    
    I = inf * uint8(ones(1500,3000,3));
    
    % stitch I_1
    indR = 30;
    indC = 810;
    I(indR:indR + size(I_1,1) - 1, indC:indC + size(I_1,2) - 1, :) = I_1;
    
    % stitch I_2
    indR = 710;
    indC = 50;
    I(indR:indR + size(I_2,1) - 1, indC:indC + size(I_2,2) - 1, :) = I_2;
    
    % stitch I_3
    indR = 710;
    indC = 950;
    I(indR:indR + size(I_3,1) - 1, indC:indC + size(I_3,2) - 1, :) = I_3;
    
    % stitch I_4
    indR = 710;
    indC = 1900;
    I(indR:indR + size(I_4,1) - 1, indC:indC + size(I_4,2) - 1, :) = I_4;
    
    figure;
    imshow(I);
    
    % print(fullfile(output_dir, ['print merged_iE=',num2str(iE),'.tiff']),'-dtiff');    
    imwrite(I, fullfile(output_dir, ['merged_iE=',num2str(iE),'.tiff']));
end
close all;

%%












