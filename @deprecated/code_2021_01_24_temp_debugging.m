% debug scripts, 2021-01-23, for twin variant identification
%% For ref. After code, the raw ID#s of the same grain at different iEs were already changed to the same
link_table = load('E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis\geotrans_and_id_link.mat', 'tbl');
link_table = link_table.tbl;

%% data for iE=0
d = load(fullfile(save_dir, ['Mg4Al_C3_parent_grain_file_iE_0.mat']));
gID_0 = d.gID;
gPhi1_0 = d.gPhi1;
gPhi_0 = d.gPhi;
gPhi2_0 = d.gPhi2;
ID_0 = d.ID;
boundary_0 = find_boundary_from_ID_matrix(ID_0);
myplot(x,y,ID_0,boundary_0);

%% (1) This is OK. Overlapped grain has too big misorientation compared to parent orientation, so not considered as twin.
iE = 1;
ii = 60;    % parent grain of interest
ID_parent = 58;
jj = 2;
ID_twin = 163;

%% (2)
iE = 2;
ii = 78;
ID_parent = 73;


%% label parent grain on map
gb_p = find_one_boundary_from_ID_matrix(ID_p);
fp = myplot(x,y,ID_p,gb_p);
label_map_with_ID(x,y,ID_p,gcf,ID_parent, 'r');

%% plot hcp_cell for parent

% id_p = gID_p(ii);
id_p = 73
ind_p = find(gID_p == id_p);
euler_p = [gPhi1_p(ind_p), gPhi_p(ind_p), gPhi2_p(ind_p)]
hcp_cell('euler', euler_p, 'ss', 25);


%% label twin grain on map
gb_c = find_one_boundary_from_ID_matrix(ID_c);
fc = myplot(x,y,ID_c,gb_c);
label_map_with_ID(x,y,ID_c,gcf, id_c, 'k', 12, 1);
% label_map_with_ID(x,y,ID_c,gcf, ID_twin, 'r', 12, 1);

%% plot hcp_cell for twin
id_twin = id_c(jj);  % twin grains id
ind_twin = (gID_c == id_twin);
euler_twin = [gPhi1_c(ind_twin), gPhi_c(ind_twin), gPhi2_c(ind_twin)];
hcp_cell('euler', euler_twin);


%%  check the parent orientation at iE = 0
% find the ID at iE=0
ind = find(gID_0==ID_parent);
euler = [gPhi1_0(ind), gPhi_0(ind), gPhi2_0(ind)]
hcp_cell('euler', euler, 'ss', 25);
calculate_misorientation_hcp(euler, euler_p)

%% check the raw orientation 
f1 = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\Mg4Al_C3 grain_file_type_1 iE=0.txt';
f2 = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\Mg4Al_C3 grain_file_type_2 iE=0.txt';
%%
f1 = ['E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\Mg4Al_C3 parent_grain_file_type_1 iE=2.txt'];
f2 = ['E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\Mg4Al_C3 parent_grain_file_type_2 iE=2.txt'];
%%
data = grain_file_to_data(f1,f2);
gb = find_one_boundary_from_ID_matrix(data.ID);
myplot(data.x, data.y, data.ID, gb)

ind = find(data.gID==73);
euler_0 = [data.gPhi1(ind), data.gPhi(ind), data.gPhi2(ind)]
hcp_cell('euler', euler_0, 'ss',25);

%% check modified orientation
f = load('E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis\Mg4Al_C3_modified_parent_grain_file_iE_2.mat');
id_in_this_step = 78;

gID = f.gID; 
gPhi1 = f.gPhi1;
gPhi = f.gPhi;
gPhi2 = f.gPhi2;

ind = find(gID==id_in_this_step)
euler = [gPhi1(ind), gPhi(ind), gPhi2(ind)]
hcp_cell('euler', euler_0, 'ss',25);




%%
iR = 420;
iC = 360;

ID_0(iR,iC)
ind_0 = find(gID_0 == 111);
euler_0 = [gPhi1_0(ind_0), gPhi_0(ind_0), gPhi2_0(ind_0)]
hcp_cell('euler', euler_0)

ID(iR,iC)
ind = find(gID == ID(iR,iC))
euler_1 = [gPhi1(ind), gPhi(ind), gPhi2(ind)]
hcp_cell('euler', euler_1)



%%
load(fullfile(save_dir, ['Mg4Al_C3_grain_file_iE_0.mat']))
ID(iR,iC)
ind = find(gID == 111)
euler_1 = [gPhi1(ind), gPhi(ind), gPhi2(ind)]
hcp_cell('euler', euler_1)

load(fullfile(save_dir, ['Mg4Al_C3_grain_file_iE_1.mat']))
ID(iR,iC)
ind = find(gID == 306)
euler_1 = [gPhi1(ind), gPhi(ind), gPhi2(ind)]
hcp_cell('euler', euler_1)

load(fullfile(save_dir, ['Mg4Al_C3_modified_parent_grain_file_iE_1.mat']))
ID(iR,iC)
ind = find(gID == 109)
euler_1 = [gPhi1(ind), gPhi(ind), gPhi2(ind)]
hcp_cell('euler', euler_1)


%%
gb = find_one_boundary_from_ID_matrix(ID_p(1:3:end,1:3:end));
a = ID_variant(1:3:end,1:3:end);
b = ID_miso(1:3:end,1:3:end);
a(b>10) = 0;
myplotm(a, gb);
caxis([-1 7]);
