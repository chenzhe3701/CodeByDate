% 2021-09-23 recover files.
% Currently, one of the problems is that 'Recovered everything in D ...' ->
% 'B recovered useful files in D', but 'Recovered partly D ...' has some more files.  
%
% First, it is better to check if there are more files in 'Recovered partly D'.  
% Then, move all useful files to 'B recovered useful ...'
% Then, move 'B recovered useful ...' useful files to 'D'.

target_dir = 'J:\Recovered everything in D as Partition [Do not moify]';
input_dir = 'J:\Recovered partly D files [Do not modify]';
output_dir = 'J:\C More files';

compare_file(input_dir, target_dir, output_dir);


function [] = compare_file(input_dir, target_dir, output_dir)

% (1) if empty, do nothing
% (2) if has folder, recursively check
% (3) if has file, check

disp(input_dir);
d = dir(input_dir);

% (1) remove '.' and '..' subfolders
ii = 1;
while(ii<=length(d))
    if strcmpi(d(ii).name,'.') || strcmpi(d(ii).name,'..')
        d(ii) = [];
    else
        ii = ii + 1;
    end
end

% recursively check
ii = 1;
while(ii<=length(d))
    input_to_check = fullfile(d(ii).folder, d(ii).name);
    target_to_check = fullfile(target_dir, d(ii).name);  % the target dir to check SHOULD have this name
    output_to_check = fullfile(output_dir, d(ii).name);
    
    % if isfolder, else isfile
    if isfolder(input_to_check)
        compare_file(input_to_check, target_to_check, output_to_check);
    else
        % if has target file, compare size; else directly copy
        if isfile(target_to_check)
            % compare file size
            if  dir(input_to_check).bytes > dir(target_to_check).bytes
                copyfile(input_to_check, output_to_check);
            end
        else
            inds = strfind(output_to_check,'\');
            ind = inds(end);
            mkdir(output_to_check(1:ind));
            copyfile(input_to_check, output_to_check);
        end
    end
    ii = ii + 1;
end

end