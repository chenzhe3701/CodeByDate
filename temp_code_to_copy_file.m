dir_name_old = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD\old grain file';
dir_name_new = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD\old grain file\New folder';
sample_name_old = 'Mg4Al_C1';
sample_name_new = 'Mg4Al_C1';
for iE = 0:7
    % iE = 0:13 for grain file
    % type 1 grain file
    old_name = fullfile(dir_name_old,[sample_name_old, '_iE=',num2str(iE),' grain_file_type_1.txt']);
    new_name = fullfile(dir_name_new,[sample_name_new, ' grain_file_type_1 iE=',num2str(iE),'.txt']);
    copyfile(old_name, new_name);
    % type 2 grain file
    old_name = fullfile(dir_name_old,[sample_name_old, '_iE=',num2str(iE),' grain_file_type_2.txt']);
    new_name = fullfile(dir_name_new,[sample_name_new, ' grain_file_type_2 iE=',num2str(iE),'.txt']);
    copyfile(old_name, new_name);
    if ismember(iE, 2:6)
        % iE = 1:13 for parent grain file
        % type 1 parent grain file
        old_name = fullfile(dir_name_old,[sample_name_old, '_iE=',num2str(iE),' grain_file_type_1_parent.txt']);
        new_name = fullfile(dir_name_new,[sample_name_new, ' parent_grain_file_type_1 iE=',num2str(iE),'.txt']);
        copyfile(old_name, new_name);
        % type 2 parent grain file
        old_name = fullfile(dir_name_old,[sample_name_old, '_iE=',num2str(iE),' grain_file_type_2_parent.txt']);
        new_name = fullfile(dir_name_new,[sample_name_new ' parent_grain_file_type_2 iE=',num2str(iE),'.txt']);
        copyfile(old_name, new_name);
    end
    
end

%% 
isfile(old_name)
