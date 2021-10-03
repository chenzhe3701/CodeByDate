% for Mohsen's Mg4Al sample, crop image, and move to corresponding folder
clc;
source_dir_p = 'E:\zhec umich Drive\2021-10-01 MgAl insitu SEM-DIC\SEM Data';

nR = 4;
nC = 5;
nT = nR * nC;   % total = nR * nC

for iS = 0:3
    source_dir = fullfile(source_dir_p, 'SE raw', ['s',num2str(iS)]);
    target_dir_SE = fullfile(source_dir_p, 'SE', ['s',num2str(iS)]);
    mkdir(target_dir_SE);
        
    for iR = 0:nR-1
        for iC = 0:nC-1
            target_dir_RC = fullfile(source_dir_p, 'RC', ['r',num2str(iR),'c',num2str(iC)]);
            mkdir(target_dir_RC);
    
            target_name = sprintf('Mg4Al_s%d_r%dc%d.tif',iS,iR,iC);
            if rem(iR,2)==0
                source_name = sprintf('image_%.5d.tif',iR*nC + iC + 1 + nT*iS);
            else
                source_name = sprintf('image_%.5d.tif',iR*nC + nC-1 -iC + 1 + nT*iS);
            end
            
            I = imread(fullfile(source_dir, source_name));
            I = imcrop(I, [0, 0, 4096, 4096]);
            
            imwrite(I, fullfile(target_dir_SE, target_name));
            imwrite(I, fullfile(target_dir_RC, target_name));
        end
    end
    
end