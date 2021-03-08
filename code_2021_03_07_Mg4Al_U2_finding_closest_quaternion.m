% Before calculating the average orientation using quaternions, we need to
% find the closest quaternion among crystallographically equivalent ones.
% This is an example discussing how to find the closest orientation.
% The example is from Mg4Al_U2, at iE=3, ii=80, parent grain ID=85, at
% step-5 (after step-4) twin analysis.
% So the deformed parent grain ID is also 85.
% child grains have IDs [421,422,424,462], where [421,462] are parent area,
% [422,424] are twinned area.

%% gather data
addChenFunction;
working_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD';
save_dir = [working_dir, '\analysis'];
mkdir(save_dir);
sample_name = 'Mg4Al_U2';

%% load data
iE= 3;
ii = 80;
id_0 = 85;
id_p = 85; % already made the same in step-4

d = load(fullfile(save_dir, 'step-4', [sample_name,'_parent_grain_file_iE_0.mat']));
gID_0 = d.gID;
gPhi1_0 = d.gPhi1;
gPhi_0 = d.gPhi;
gPhi2_0 = d.gPhi2;
ID_0 = d.ID;
boundary_0 = find_boundary_from_ID_matrix(ID_0);
phi1_0 = d.phi1;
phi_0 = d.phi;
phi2_0 = d.phi2;

d = load(fullfile(save_dir, 'step-4', [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
x = d.x;
y = d.y;
gID_p = d.gID;
gPhi1_p = d.gPhi1;
gPhi_p = d.gPhi;
gPhi2_p = d.gPhi2;
ID_p = d.ID;
boundary_p = find_boundary_from_ID_matrix(ID_p);
phi1_p = d.phi1;
phi_p = d.phi;
phi2_p = d.phi2;

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

%% plot IPF map at iE=0
ind_local = ismember(ID_0, id_0);
indC_min = find(sum(ind_local, 1), 1, 'first') - 10;
indC_max = find(sum(ind_local, 1), 1, 'last') + 10;
indR_min = find(sum(ind_local, 2), 1, 'first') - 10;
indR_max = find(sum(ind_local, 2), 1, 'last') + 10;

ID_0_local = ID_0(indR_min:indR_max, indC_min:indC_max);
boundary_0_local = find_one_boundary_from_ID_matrix(ID_0_local);
phi1_0_local = phi1_0(indR_min:indR_max, indC_min:indC_max);
phi_0_local = phi_0(indR_min:indR_max, indC_min:indC_max);
phi2_0_local = phi2_0(indR_min:indR_max, indC_min:indC_max);
x_local = x(indR_min:indR_max, indC_min:indC_max);
y_local = y(indR_min:indR_max, indC_min:indC_max);
IPF = generate_IPF_map_pixel_wise(x_local, y_local, phi1_0_local, phi_0_local, phi2_0_local, boundary_0_local, [0 0 0], [1 0 0]);
figure;
image(unique(x_local(:)), unique(y_local(:)), IPF);
axis equal;
label_map_with_ID(x_local, y_local, ID_0_local, gcf, 85,'k',18,1);

%% Plot IPF map at iE=3 for parent
ind_local = ismember(ID_p, id_p);
indC_min = find(sum(ind_local, 1), 1, 'first') - 10;
indC_max = find(sum(ind_local, 1), 1, 'last') + 10;
indR_min = find(sum(ind_local, 2), 1, 'first') - 10;
indR_max = find(sum(ind_local, 2), 1, 'last') + 10;

ID_p_local = ID_p(indR_min:indR_max, indC_min:indC_max);
boundary_p_local = find_one_boundary_from_ID_matrix(ID_p_local);
phi1_p_local = phi1_p(indR_min:indR_max, indC_min:indC_max);
phi_p_local = phi_p(indR_min:indR_max, indC_min:indC_max);
phi2_p_local = phi2_p(indR_min:indR_max, indC_min:indC_max);
x_local = x(indR_min:indR_max, indC_min:indC_max);
y_local = y(indR_min:indR_max, indC_min:indC_max);
IPF = generate_IPF_map_pixel_wise(x_local, y_local, phi1_p_local, phi_p_local, phi2_p_local, boundary_p_local, [0 0 0], [1 0 0]);
figure;
image(unique(x_local(:)), unique(y_local(:)), IPF);
axis equal;
label_map_with_ID(x_local, y_local, ID_p_local, gcf, 85,'k',18,1);

%% Plot IPF map at iE=3 for children
ID_c_local = ID_c(indR_min:indR_max, indC_min:indC_max);
boundary_c_local = find_one_boundary_from_ID_matrix(ID_c_local);
phi1_c_local = phi1_c(indR_min:indR_max, indC_min:indC_max);
phi_c_local = phi_c(indR_min:indR_max, indC_min:indC_max);
phi2_c_local = phi2_c(indR_min:indR_max, indC_min:indC_max);
IPF = generate_IPF_map_pixel_wise(x_local, y_local, phi1_c_local, phi_c_local, phi2_c_local, boundary_c_local, [0 0 0], [1 0 0]);
figure;
image(unique(x_local(:)), unique(y_local(:)), IPF);
axis equal;
label_map_with_ID(x_local, y_local, ID_c_local, gcf, [421,422,424,462],'k',18,1);

%% get data for the grains, and run codes for analysis
close all;
clc;

ii = 80;
id_p = gID_p(ii)
ind_0 = find(gID_0 == id_p);
euler_0 = [gPhi1_0(ind_0), gPhi_0(ind_0), gPhi2_0(ind_0)]
hcp_cell('euler', euler_0);

id_c = unique(ID_c(ID_p == id_p))
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
tbl = table(col_ID, col_euler, misorientation)
    
inds = find(misorientation < 15);
id_po = id_c(inds)

% child grain with parent orientation
ind = find(ismember(gID_c, 421));
euler_421 = [gPhi1_c(ind), gPhi_c(ind), gPhi2_c(ind)]
ind = find(ismember(gID_c, 462));
euler_462 = [gPhi1_c(ind), gPhi_c(ind), gPhi2_c(ind)]
% calculated parent orientation
euler_po = calculate_average_dominant_euler_hcp([euler_421; euler_462])

% child grain with twin orientation
ind = find(ismember(gID_c, 422));
euler_422 = [gPhi1_c(ind), gPhi_c(ind), gPhi2_c(ind)]
misorientation_422 = [];
for kk = 1:6
    euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
    misorientation_422(kk) = calculate_misorientation_euler_d(euler_422, euler_kk, 'HCP');
end
display(misorientation_422);

ind = find(ismember(gID_c, 424));
euler_424 = [gPhi1_c(ind), gPhi_c(ind), gPhi2_c(ind)]
misorientation_424 = [];
for kk = 1:6
    euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
    misorientation_424(kk) = calculate_misorientation_euler_d(euler_424, euler_kk, 'HCP');
end
display(misorientation_424);

%% summary, the eulers we need to use.
% parent grain at iE=0
euler_0 = [13.5330  161.0100 -179.9000];
% orientation of child grains with parent orientation at iE=3
euler_421 = [3.5090  161.5360  169.9800];
euler_462 = [6.7330  161.0710  173.4700];
% orientation of child grains, which should be identified as twins (should
% have about 87 misorientation with parent orientation);
euler_422 = [-103.7350  100.8480   77.8400];
euler_424 = [75.7620   87.3340 -134.5500];


%% Find average orientation of 421 and 462.
close all;
euler_input = [euler_421; euler_462];
[~,q0,q1,q2,q3] = euler_to_quaternion(euler_input(:,1), euler_input(:,2), euler_input(:,3));
qs = [q0,q1,q2,q3];
q_ref = qs(1,:);
q_input = qs(2,:);

% (step) find closest orientation among equivalent orientation
S = [1, 0, 0, 0;
    sqrt(3)/2, 0, 0, 1/2;
    1/2, 0, 0, sqrt(3)/2;
    0, 0, 0, -1;
    1/2, 0, 0, -sqrt(3)/2;
    sqrt(3)/2, 0, 0, -1/2;
    0, 1, 0, 0;
    0, sqrt(3)/2, -1/2, 0;
    0, 1/2, -sqrt(3)/2, 0;
    0, 0, 1, 0;
    0, 1/2, sqrt(3)/2, 0;
    0, sqrt(3)/2, 1/2, 0];
% crystallographically equivalent ones
q_equivalent = quatmultiply(q_input, S);

% calculate difference
q_delta = quatmultiply(q_equivalent, quatconj(q_ref));
axis_angle = quat2axang(q_delta); 
thetad = axis_angle(:,4)/pi*180;
[val_min, ind] = min(abs(thetad));
q_closest = q_equivalent(ind,:);

method = 1; % 1 = make axis similar, 2 = make angle the same, 3 = 
% analyze q_out, see if want to take negative
if method == 1
    % make axis similar
    if dot(q_closest(2:4),q_ref(2:4)) < 0
        q_closest = -q_closest;
    end    
elseif method == 2
    % make angle same sign
    if  q_closest(1) * q_ref(1) < 0
        q_closest = -q_closest;
    end
end

% The eulers (and M) is not dependent on the sign of the quaternion.
m = quat2dcm(q_closest);
[a,b,c] = dcm2angle(m,'zxz');
euler_closest = [a,b,c];
euler_closest = euler_closest/pi*180


% (step) calculate average
qavg = calculate_avg_quat([q_ref;q_closest]);
m = quat2dcm(qavg);
[a,b,c] = dcm2angle(m,'zxz');
euler_avg = [a,b,c];
euler_avg = euler_avg/pi*180;

% plot unit cell
hcp_cell('euler', euler_421);
hcp_cell('euler', euler_462);
hcp_cell('euler', euler_avg);
