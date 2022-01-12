%% Move photos and movies, from source dir, to target dir
clear;
clc;

% source dirs
input_dirs = {'D:\D Folder\101MSDCF\iphone 2011';
    'D:\D Folder\101MSDCF\iphone2 2011-2012';
    'D:\D Folder\101MSDCF\iphone3 2012';
    'D:\D Folder\101MSDCF\iphone4 2014-2015';
    'D:\D Folder\101MSDCF\iphone5 2016-2017';
    'D:\D Folder\101MSDCF\iphone6 2017-2018';
    'D:\D Folder\101MSDCF\iphone6 2017-2018';
    'C:\Users\ZheChen\Pictures\iCloud Photos\Downloads'};

do_copy = 1;
% target dir
output_dir = fullfile('D:','photos');
mkdir(output_dir);
mkdir(fullfile(output_dir, 'Non Images'));
mkdir(fullfile(output_dir, 'Other Model'));
mkdir(fullfile(output_dir, 'iPhone'));

% for each folder: find subfolders
file_structure = [];
for ii = 1:length(input_dirs)
    file_structure = [file_structure; dir(fullfile(input_dirs{ii},'/**/*.*'))];
end

% remove directories from the structure
for ii = length(file_structure):-1:1
    if file_structure(ii).isdir
        file_structure(ii) = [];
    end
end

% remove non-images from 'file_structure', copy non-images to target folder
for ii = length(file_structure):-1:1
    try 
        % read imf info and append
        img_info = imfinfo(fullfile(file_structure(ii).folder, file_structure(ii).name));
        img_info_field_names = fieldnames(img_info);
        
        % copy field from [img_info] to [file_structure]
        for jj = 1:length(img_info_field_names)
            img_info_field_name = img_info_field_names{jj};
            file_structure(ii).(img_info_field_name) = img_info.(img_info_field_name);
        end
    
    catch
        % remove non-images from [file_structure], and copy to non-image folder  
        if do_copy
            copyfile(fullfile(file_structure(ii).folder, file_structure(ii).name), fullfile(output_dir, 'Non Images', file_structure(ii).name));
        end
        file_structure(ii) = [];
    end
end

% 1st sort, by datenum
ind_1 = 1;
count_redundant_1 = 0;
N = size(file_structure,1);
indr_ok = zeros(N,1);   % unique rows
while ind_1 <= N
   % (1) check redundancy, move ind_1 to the last of the same image 
    ind_2 = ind_1;
   % while row 1 and row 2 have same datenum
   while ind_2 < N && (file_structure(ind_1).datenum == file_structure(ind_2+1).datenum) && (file_structure(ind_1).bytes ==file_structure(ind_2+1).bytes)
      ind_2 = ind_2 + 1;
      count_redundant_1 = count_redundant_1 + 1;
   end
   
   if ind_2 > ind_1
       ind_1 = ind_2;
   end
   indr_ok(ind_1) = 1;
   ind_1 = ind_1 + 1;
end
file_structure = file_structure(logical(indr_ok),:);

% 2nd sort, by 'DateTime'
tbl = struct2table(file_structure);
for i = 1:size(tbl,1)
    try
        tbl.SortCol(i) = datenum(tbl.DateTime{i}, 'yyyy:mm:dd HH:MM:SS');
    catch
        tbl.SortCol(i) = tbl.datenum(i);
    end
end
[tbl,idx] = sortrows(tbl, {'SortCol','name'});
file_structure = table2struct(tbl);    % file_structure(idx);

ind_1 = 1;
count_redundant_2 = 0;
N = size(file_structure,1);
indr_ok = zeros(N,1);   % unique rows
while ind_1 <= N
   % (1) check redundancy, move ind_1 to the last of the same image 
    ind_2 = ind_1;
   % while row 1 and row 2 have same datenum
   while ind_2 < N && (file_structure(ind_1).SortCol == file_structure(ind_2+1).SortCol) && (file_structure(ind_1).bytes ==file_structure(ind_2+1).bytes)
      ind_2 = ind_2 + 1;
      count_redundant_2 = count_redundant_2 + 1;
   end
   
   if ind_2 > ind_1
       ind_1 = ind_2;
   end
   indr_ok(ind_1) = 1;
   ind_1 = ind_1 + 1;
end
file_structure = file_structure(logical(indr_ok),:);


% For the remaining images, copy to targat folder after checking redundant, target based on model name   
ind_1 = 1;
count_A = 1;
count_B = 1;

N = size(file_structure,1);
while ind_1 <= N

   % determine [target_dir] and [target_name] based on [mdl] and [fname]
   model_name = file_structure(ind_1).Model;
   file_name = file_structure(ind_1).name;
   if ~isempty(model_name) && contains(model_name, 'iPhone')
       % determine [target_dir]
       target_dir = fullfile(output_dir, 'iPhone');   
       % determine [target_name]
       if regexp(file_name, 'IMG_')
           [~,ext] = strtok(file_name,'.');
           target_name = ['IMG_',num2str(count_A, '%.4d'),ext];
           count_A = count_A + 1;
       else
           target_name = file_name;
           % if name exist, append a letter, such as '_A', '_B', etc
           iLetter = 65;
           while exist(fullfile(target_dir, target_name), 'file')
               [a,b] = strtok(target_name,'.');
               target_name = [a,'_',char(iLetter),b];
               iLetter = iLetter + 1;
           end
       end
   else
       % determine [target_dir]
       target_dir = fullfile(output_dir, 'Other Model');
       % determine [target_name]
       if regexp(file_name, 'IMG_')
           [~,ext] = strtok(file_name,'.');
           target_name = ['IMG_',num2str(count_B, '%.4d'),ext];
           count_B = count_B + 1;
       else
           target_name = file_name;
           % if name exist, append a letter, such as '_A', '_B', etc
           iLetter = 65;
           while exist(fullfile(target_dir, target_name), 'file')
               [a,b] = strtok(target_name,'.');
               target_name = [a,'_',char(iLetter),b];
               iLetter = iLetter + 1;
           end
       end
   end
   
   % copy file
   if do_copy
       copyfile(fullfile(file_structure(ind_1).folder, file_name), fullfile(target_dir, target_name));
   end
   % then increment ind_1
   ind_1 = ind_1 + 1;
end

   
   
  
    
    
    
    
    
    
    
    
    
    
    
    


