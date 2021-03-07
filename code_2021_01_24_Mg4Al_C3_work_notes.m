% chenzhe, 2021-01-24
% problem: It seems that for Mg4Al_C3, the raw EBSD identified grain
% orientation is different.
% The reason might be:
% If at iE=2, the dominant part of this grain is already twinned, then the
% dominant orientation is the 'twin orientation' rather than 'parent
% orientation', because OIM does not know the parent information at iE=0.
%
% This should have been taken care of in my updated twin variant
% identification algorithm, i.e., (1) check parent orienation at iE=0, (2)
% check children with parent orientation at iE>0, take their average, and
% use that to represent the parent orientation at iE>0.
%% at iE=0, grain id = 73
f1 = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\Mg4Al_C3 grain_file_type_1 iE=0.txt';
f2 = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\Mg4Al_C3 grain_file_type_2 iE=0.txt';

data = grain_file_to_data(f1,f2);
gb = find_one_boundary_from_ID_matrix(data.ID);
myplot(data.x, data.y, data.ID, gb);
label_map_with_ID(data.x, data.y, data.ID, gcf, 73, 'r');

ind = find(data.gID==73);
euler_0 = [data.gPhi1(ind), data.gPhi(ind), data.gPhi2(ind)]
hcp_cell('euler', euler_0, 'ss',25);
title('raw euler before aligning to sample');

%% at iE=0, grain id = 78
f1 = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\Mg4Al_C3 parent_grain_file_type_1 iE=2.txt';
f2 = 'E:\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\Mg4Al_C3 parent_grain_file_type_2 iE=2.txt';

data = grain_file_to_data(f1,f2);
gb = find_one_boundary_from_ID_matrix(data.ID);
myplot(data.x, data.y, data.ID, gb);
label_map_with_ID(data.x, data.y, data.ID, gcf, 78, 'r');

ind = find(data.gID==78);
euler = [data.gPhi1(ind), data.gPhi(ind), data.gPhi2(ind)]
hcp_cell('euler', euler, 'ss',25);
title('raw euler before aligning to sample');