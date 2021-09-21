clear;clc;

% copy and rename images into appropriate folders
main_folder = 'E:\zhec umich Drive\2021-09-11 UM134_Mg_D3 insitu SEM-DIC';
sample_name = 'UM134_Mg_D3';

%% built-in images
n_0 = [1,1,37,73,109,145,181,1];    % the first number of built-in images taken at that iE

nR = 6;
nC = 6;

for iE = 0:7
    iB = iE + 1;
    
    source_dir = fullfile(main_folder, 'raw', ['iE=',num2str(iE)]);
    target_dir = fullfile(main_folder, 'SEM Data', 'SE builtin', ['s',num2str(iE)]);
    mkdir(target_dir);       

    for iR = 0:(nR-1)
        for iC = 0:(nC-1)
            if(rem(iR,2)==0)
                k = iR*nC + iC + n_0(iB);
            else
                k = iR*nC + (nC-1-iC) + n_0(iB);
            end

            source_file = ['image_',num2str(k,'%05d'),'.tif'];
            target_file = [sample_name, '_s',num2str(iE), '_r',num2str(iR),'c',num2str(iC),'.tif'];
            
            I = imread(fullfile(source_dir,source_file));
            I = imcrop(I, [0 0 4096 4096]);
            imwrite(I, fullfile(target_dir,target_file));
            
        end
    end
    
end

%% external images
ext_name_part = {'','_s1','_s2','_s3','_s4','_s5','_6','_s7'};  % external images names are not consistent ... 
nR = 21;
nC = 21;

for iE = 0:7
    iB = iE + 1;
    
    source_dir = fullfile(main_folder, 'raw', ['iE=',num2str(iE)]);
    target_dir = fullfile(main_folder, 'SEM Data', 'SE external', ['s',num2str(iE)]);
    mkdir(target_dir);       

    for iR = 0:(nR-1)
        for iC = 0:(nC-1)           
            source_file = [sample_name, ext_name_part{iB}, '_r',num2str(iR),'c',num2str(iC),'.tiff'];
            target_file = [sample_name, '_s',num2str(iE), '_r',num2str(iR),'c',num2str(iC),'.tif'];
            copyfile(fullfile(source_dir,source_file), fullfile(target_dir,target_file), 'f');
        end
    end
    
end