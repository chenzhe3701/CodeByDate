% study the meaning of 'anti-grains' in OIM data.
% The data was from Mg4Al_C1 tested on 2020-10-23, iE=0
% Data was Neighbor Orientation Correlation (NOC) cleaned with tolerange
% angle = 5, minimum CI = 0.1, clean up level = 3.
% ------
% Conclusion: 
% (1) Grains not in partition will be numbered as '0'.
% (2) Anti-grains in partition have IDs, which are > than IDs of valid grains. 

%% load data and show header
addChenFunction;
cdi(6);
[EBSD_data_1, EBSD_header_1] = grain_file_read('data\for_anti_grain_file_type_1.txt');
[EBSD_data_2, EBSD_header_2] = grain_file_read('data\for_anti_grain_file_type_2.txt');
%% find column
column_index_1 = find_variable_column_from_grain_file_header(EBSD_header_1,...
    {'grain-ID','x-um','y-um'});
column_index_2 = find_variable_column_from_grain_file_header(EBSD_header_2,...
    {'grainId'});

%% read type-2 grain file and get average info for grains
gID = EBSD_data_2(:,column_index_2(1));

%% EBSD data, from type-1 grain file. (column, data) pair:
% (1,phi1) (2,phi) (3,phi2) (4,xMicron) (5,yMicron) (6,IQ) (7,CI) (8,Fit) (9,grain-Id) (10,edgeGrain?)
% Read EBSD data.  IQ,CI,Fit are not needed for now, but might need in future
x = EBSD_data_1(:,column_index_1(2));
y = EBSD_data_1(:,column_index_1(3));
unique_x = unique(x(:));
ebsdStepSize = unique_x(2) - unique_x(1);
mResize = (max(x(:)) - min(x(:)))/ebsdStepSize + 1;
nResize = (max(y(:)) - min(y(:)))/ebsdStepSize + 1;

ID = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
%%
min(gID(:))
max(gID(:))
min(ID(:))
max(ID(:))
myplot(ID);
