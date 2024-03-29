% try to stitch images for SEM-DIC sample UM134_Mg_D3

% DRAWFIGURE = 0: no figure, 1: some, 2: more
DRAWFIGURE = 0;

% Path for images and dic data (dic data usually in the same folder as the images)
path_DIC = 'E:\zhec umich Drive\2021-09-11 UM134_Mg_D3 insitu SEM-DIC\SEM Data\SE external';
path_target = 'E:\zhec umich Drive\2021-09-11 UM134_Mg_D3 insitu SEM-DIC\SEM Data\stitched images external';

% subfolder name, 's0', 's1', ...
subfolderNamePrefix_1 = 's';

fileNamePrefix_1 = 'UM134_Mg_D3_s';  
fileNamePrefix_2 = '_';

% resolution of images
resX = 1024;
resY = 1024;
% overlay/window size to search and match images
OVERLAY = 200;
% save reduced size image
reduction = 10;

B = 1;   % 'B' for 'base', to handle if it's 0/1-based index.  But B=1 for 0-based. B=0 for 1-based.  When iR, iC is used with FOV, transX, ... add this B.
row_start = 0; % starting # of FOV rows
row_end = 20;
col_start = 0;
col_end = 20;  % ending # of FOV cols
iE_start = 0;
iE_stop = 7;
singleFOV = 0;  % This overwrite the 'iE' so the code can be conveniently applied to a single strain
singleRow = 0;

% FOV = 'r0c0', ...
FOV = make_FOV_string(abs(B-1), row_end, abs(B-1), col_end, 1, 'rc');   

% initialize
clear transX; clear transY; clear transX_incremental; clear transY_incremental;
transX_incremental = zeros(row_end+B,col_end+B);
transY_incremental = zeros(row_end+B,col_end+B);
transX = zeros(row_end+B,col_end+B);
transY = zeros(row_end+B,col_end+B);
clear specialRC;    % but can define special cases
corrMethod = 3;     % 1 = fft, 2 = normxcorr2
cutEdge = 0;    % cut edge = 1, vs average=0, in blending


%%
for iE = iE_start:iE_stop        % 'e#' in the file name, i.e., stop/pause #  ----------------------------------------------
    if singleFOV == 1
       iE = []; 
    end
    subfolderName = [subfolderNamePrefix_1,num2str(iE)];
    Oly = OVERLAY;
    for iC = col_start:col_end
        for iR = row_start:row_end - 1 + singleRow
            % search the down-sided picture
            pause(1);
            close all;
            try
                fName1 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,FOV{iR+B,iC+B}]; disp(fName1);
                fName2 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,FOV{iR+B+1,iC+B}]; disp(fName2);
                I = imread([path_DIC,'\',subfolderName,'\',fName1,'.tif']);
                J = imread([path_DIC,'\',subfolderName,'\',fName2,'.tif']);
                stitch_direction = 1;   % 1 for up-down, 2 for left-right
                % For debug
                if DRAWFIGURE > 1
                    f1 = figure;imshowpair(I,J,'montage');
                end
                
                try
                    switch corrMethod
                        case 1
                            [yOffSet,xOffSet] = fft_register(I,J,'d',[0.7*size(I,1)*0, 0,  Oly*0, 0], [0, 0.7*size(J,1)*0,  Oly*0, 0], 1);         % could change this accordingly, be careful with your choice of parameter ---------------------------
                        case 2
                            % initially crop a small region to detect
                            yi = size(I,1) - Oly; % xi, yi are zero-based.
                            xi = round(size(I,2)/2-Oly/2);
                            
                            Iprime = imcrop(I,[xi,yi,Oly,Oly]);
                            x_neg = size(I,2)-xi;   % Iprime wrt lower right corner of I
                            y_neg = size(I,1)-yi;
                            
                            c = fft_xcorr2(double(J), double(Iprime), 1);
                            
                            if DRAWFIGURE > 1
                                figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
                            end
                            [ypeak, xpeak] = find(c == max(c(:)));
                            xOffSet_J = xpeak-size(Iprime,2);
                            yOffSet_J = ypeak-size(Iprime,1);
                            
                            osx_neg = x_neg + xOffSet_J;    % OffSet_Negative, J's upper left corner's offset wrt I's lower-right corner
                            osy_neg = y_neg + yOffSet_J;
                            
                            xOffSet = size(I,2) - osx_neg;  % positive offset, J's upper-left coner's offset wrt I's upper-left corner
                            yOffSet = size(I,1) - osy_neg;
                        case 3
                            [yOffSet,xOffSet] = normxcorr2A_register(J, I, [0, 0.7*size(J,1),  Oly*0.1, Oly*0.1], [0.7*size(I,1), 0,  Oly*0.1, Oly*0.1], 1)
                    end
                    
                    if DRAWFIGURE > 0
                        figure;imshowpair(I,J,'montage');
                        imrect(gca, [xi, yi, size(Iprime,2), size(Iprime,1)]);
                        imrect(gca, [size(I,2)+xOffSet_J, yOffSet_J, size(Iprime,2), size(Iprime,1)]);
                    end
                    
                    if (stitch_direction==1) && ( (yOffSet > size(I,1))||(abs(xOffSet)>500) )
                        error();
                    end
                    
                catch
                    % get rectangular region manually
                    try close(f1); catch ; end
                    f1 = figure;set(f1,'name','select an area to match, prefer on image I (left)');imshowpair(I,J,'montage');
                    rect = round(getrect(gcf));
                    xi = rect(1); yi = rect(2); OlyX = rect(3); OlyY = rect(4);
                    xi = mod(xi-size(I,2),size(I,2));
                    OlyX = mod(OlyX-size(I,2),size(I,2));
                    yi = mod(yi-size(I,1),size(I,1));
                    OlyY = mod(OlyY-size(I,1),size(I,1));
                    Iprime = imcrop(I,[xi,yi,OlyX,OlyY]);
                    x_neg = size(I,2)-xi;   % Iprime wrt lower right corner of I
                    y_neg = size(I,1)-yi;
                    
                    c = fft_xcorr2(double(J), double(Iprime), 1);
                    if DRAWFIGURE > 1
                        figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
                    end
                    [ypeak, xpeak] = find(c == max(c(:)));
                    xOffSet_J = xpeak-size(Iprime,2);
                    yOffSet_J = ypeak-size(Iprime,1);
                    
                    osx_neg = x_neg + xOffSet_J;    % J's upper left corner's offset wrt I's lower-right corner
                    osy_neg = y_neg + yOffSet_J;
                    
                    xOffSet = size(I,2) - osx_neg;
                    yOffSet = size(I,1) - osy_neg;
                    
                    if DRAWFIGURE > 0
                        figure;imshowpair(I,J,'montage');
                        imrect(gca, [xi, yi, size(Iprime,2), size(Iprime,1)]);
                        imrect(gca, [size(I,2)+xOffSet_J, yOffSet_J, size(Iprime,2), size(Iprime,1)]);
                    end
                    
                end
                
                transX_incremental(iR+B+1,iC+B) = xOffSet
                transY_incremental(iR+B+1,iC+B) = yOffSet
                transX(iR+B+1,iC+B) = xOffSet  + transX(iR+B,iC+B);
                transY(iR+B+1,iC+B) = yOffSet  + transY(iR+B,iC+B);
            end
            
            if (iR==row_start)&&(iC<col_end)    % search the right-side picture as well
                pause(1);
                close all;
                fName1 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,FOV{iR+B,iC+B}]; disp(fName1);
                fName3 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,FOV{iR+B,iC+B+1}]; disp(fName3);
                I = imread([path_DIC,'\',subfolderName,'\',fName1,'.tif']);
                J = imread([path_DIC,'\',subfolderName,'\',fName3,'.tif']);
                stitch_direction = 2;   % 1 for up-down, 2 for left-right
                % For debug
                if DRAWFIGURE > 1
                    figure;imshowpair(I,J,'montage');
                end
                
                try
                    switch corrMethod
                        case 1
                            [yOffSet,xOffSet] = fft_register(I,J,[Oly*0, 0,  0.7*size(I,2)*0, 0], [0, Oly*0,  0, 0.7*size(J,2)*0], 1);       % change this accordingly, be careful with your choice of parameter ---------------------------
                        case 2
                            % initially crop a small region to detect
                            xi = size(I,2)-Oly; % xi, yi are zero-based.
                            
                            if iR == row_start
                                yi = size(I,1)-Oly;     % This is a value that can be fine tuned. E.g., add '-300' to treat Rotated data for test '20170430' -----------------
                            elseif iR == row_end
                                yi = 0;                 % This is a value that can be fine-tuned. E.g., add '+100' to treat Rotated data for test '20170430'------------------
                            else
                                yi = round(size(I,1)/2-Oly/2);
                            end
                            
                            Iprime = imcrop(I,[xi,yi,Oly,Oly]);
                            x_neg = size(I,2)-xi;   % Iprime wrt lower right corner of I
                            y_neg = size(I,1)-yi;
                            
                            c = fft_xcorr2(double(J), double(Iprime), 1);
                            
                            if DRAWFIGURE > 1
                                figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
                            end
                            [ypeak, xpeak] = find(c == max(c(:)));
                            xOffSet_J = xpeak-size(Iprime,2);
                            yOffSet_J = ypeak-size(Iprime,1);
                            
                            osx_neg = x_neg + xOffSet_J;    % J's upper left corner's offset wrt I's lower-right corner
                            osy_neg = y_neg + yOffSet_J;
                            
                            xOffSet = size(I,2) - osx_neg;
                            yOffSet = size(I,1) - osy_neg;
                        case 3
                            [yOffSet,xOffSet] = normxcorr2A_register(J, I, [Oly*0.1, Oly*0.1,  0, 0.7*size(J,2)], [Oly*0.1, Oly*0.1,  0.7*size(I,2), 0], 1);
                            
                    end
                    
                    if DRAWFIGURE > 0
                        figure;imshowpair(I,J,'montage');
                        imrect(gca, [xi, yi, size(Iprime,2), size(Iprime,1)]);
                        imrect(gca, [size(I,2)+xOffSet_J, yOffSet_J, size(Iprime,2), size(Iprime,1)]);
                     end
                    
                    if (stitch_direction==2) && ( (xOffSet > size(I,2))||(abs(yOffSet)>500) )
                        error();
                    end
                    
                catch
                    % get rectangular region manually
                    try close(f1); catch ; end
                    f1 = figure;set(f1,'name','select an area to match, prefer on image I (left)');imshowpair(I,J,'montage');
                    rect = round(getrect(gcf));
                    xi = rect(1); yi = rect(2); OlyX = rect(3); OlyY = rect(4);
                    xi = mod(xi-size(I,2),size(I,2));
                    OlyX = mod(OlyX-size(I,2),size(I,2));
                    yi = mod(yi-size(I,1),size(I,1));
                    OlyY = mod(OlyY-size(I,1),size(I,1));
                    Iprime = imcrop(I,[xi,yi,OlyX,OlyY]);
                    x_neg = size(I,2)-xi;   % Iprime wrt lower right corner of I
                    y_neg = size(I,1)-yi;
                    
                    c = fft_xcorr2(double(J), double(Iprime), 1);
                    if DRAWFIGURE > 1
                        figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
                    end
                    [ypeak, xpeak] = find(c == max(c(:)));
                    xOffSet_J = xpeak-size(Iprime,2);
                    yOffSet_J = ypeak-size(Iprime,1);
                    
                    osx_neg = x_neg + xOffSet_J;    % J's upper left corner's offset wrt I's lower-right corner
                    osy_neg = y_neg + yOffSet_J;
                    
                    xOffSet = size(I,2) - osx_neg;
                    yOffSet = size(I,1) - osy_neg;
                    
                    if DRAWFIGURE > 0
                        figure;imshowpair(I,J,'montage');
                        imrect(gca, [xi, yi, size(Iprime,2), size(Iprime,1)]);
                        imrect(gca, [size(I,2)+xOffSet_J, yOffSet_J, size(Iprime,2), size(Iprime,1)]);
                    end
                    
                end
                
                transX_incremental(iR+B,iC+B+1) = xOffSet
                transY_incremental(iR+B,iC+B+1) = yOffSet
                transX(iR+B,iC+B+1) = xOffSet + transX(iR+B,iC+B);
                transY(iR+B,iC+B+1) = yOffSet + transY(iR+B,iC+B);
                
            end
            
            % if have specialRC to handle [currentR,currentC, neighborR,neighborC, assignedTransX, assignedTransY]
            ind = [];
            if exist('specialRC','var')&&(~isempty(specialRC))
                ind = find((iR==specialRC(:,1))&(iC==specialRC(:,2))&(iR==specialRC(:,3))&(iC+1==specialRC(:,4)));  % right neighbor
                if isempty(ind)
                    ind = find((iR==specialRC(:,1))&(iC==specialRC(:,2))&(iR+1==specialRC(:,3))&(iC==specialRC(:,4)));  % bottom neighbor
                end
            end
            if ~isempty(ind)
                transX_incremental(specialRC(ind,3)+B, specialRC(ind,4)+B) = specialRC(ind,5)
                transY_incremental(specialRC(ind,3)+B, specialRC(ind,4)+B) = specialRC(ind,6)
                transX(specialRC(ind,3)+B, specialRC(ind,4)+B) = specialRC(ind,5) + transX(specialRC(ind,1)+B,specialRC(ind,2)+B);
                transY(specialRC(ind,3)+B, specialRC(ind,4)+B) = specialRC(ind,6) + transY(specialRC(ind,1)+B,specialRC(ind,2)+B);
            end
            
        end
    end
    
    makeAllTransPositive = 1;
    if makeAllTransPositive
        xi = -min(transX(:));     % overall shift when upper_left is (0,0) --------------------------
        yi = -min(transY(:));
        transX = transX + xi;   % make the smallest translation as 0, instead of negative
        transY = transY + yi;
        transX_incremental(row_start+B,col_start+B) = transX_incremental(row_start+B,col_start+B) + xi;
        transY_incremental(row_start+B,col_start+B) = transY_incremental(row_start+B,col_start+B) + yi;
    end
    
    save([path_target,'\','translations_searched_vertical_stop_',num2str(iE)],'transX','transY','transX_incremental','transY_incremental');
    
    % The following is the stitching part.
    xi = -min(transX(:));     % overall shift when upper_left is (0,0) --------------------------
    yi = -min(transY(:));
    xm = max(transX(:));    % maximum shift
    ym = max(transY(:));
    img = zeros(resY+yi+ym,resX+xi+xm) * nan;
    
    for iR = row_start:row_end
        for iC = col_start:col_end
            fName = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,FOV{iR+B,iC+B}]; disp(fName);    % change this accordingly ----------------------------------------------------------
            
            I = double(imread([path_DIC,'\',subfolderName,'\',fName,'.tif']));    % otherwise, nanmean doesn't work, bucause cat makes nan to 0.
            
            % cutted row/col index = 1: (overlap_boundary-overlay_size)
            if cutEdge
                if (iR>row_start)&&(iC>col_start)
                    I(1:size(I,1)-(transY(iR+B,iC+B)-transY(iR+B-1,iC+B))-Oly, :) = nan;    % cut top, also consider top neighbor
                    I(:, 1:size(I,2)-(transX(iR+B,iC+B)-transX(iR+B,iC+B-1))-Oly) = nan;    % cut left, also consider left neighbor
                elseif (iC>col_start)
                    I(:, 1:size(I,2)-(transX(iR+B,iC+B)-transX(iR+B,iC+B-1))-Oly) = nan;    % iR = first row, cut left, only consider left neighbor
                elseif (iR>row_start)
                    I(1:size(I,1)-(transY(iR+B,iC+B)-transY(iR+B-1,iC+B))-Oly, :) = nan;    % iC = first col, cut top, only consider top neighbor
                end
            end
            
            J = img(1+transY(iR+B,iC+B)+yi:resY+transY(iR+B,iC+B)+yi,1+transX(iR+B,iC+B)+xi:resX+transX(iR+B,iC+B)+xi);
            
            img(1+transY(iR+B,iC+B)+yi:resY+transY(iR+B,iC+B)+yi,1+transX(iR+B,iC+B)+xi:resX+transX(iR+B,iC+B)+xi) = mean(cat(3,I,J),3,'omitnan');
        end
    end
    
    bs = whos('img');
    bs = bs.bytes;
    if bs > (2^32-1)*4
        img = uint8(img/255);
        imwrite(img,[path_target,'\img_full_res_iE=',num2str(iE),'.tif'],'tif');
        imwrite(img(1:2:end,1:2:end),[path_target,'\img_half_res_iE=',num2str(iE),'.tif'],'tif');
        imwrite(img(1:10:end,1:10:end),[path_target,'\img_10_res_iE=',num2str(iE),'.tif'],'tif');

    else
        img = uint16(img);
        imwrite(img,[path_target,'\img_full_res_iE=',num2str(iE),'.tif'],'tif');
        imwrite(img(1:2:end,1:2:end),[path_target,'\img_half_res_iE=',num2str(iE),'.tif'],'tif');
        imwrite(img(1:10:end,1:10:end),[path_target,'\img_10_res_iE=',num2str(iE),'.tif'],'tif');
        
%         % a code that can divide image into 4 smaller images. disabled but kept code.
%         imwrite(img(1:end/2,1:end/2),[path_target,'\fullImage_1_',num2str(iE),'.tif'],'tif');
%         imwrite(img(1:end/2,end/2+1:end),[path_target,'\fullImage_2_',num2str(iE),'.tif'],'tif');
%         imwrite(img(end/2+1:end,1:end/2),[path_target,'\fullImage_3_',num2str(iE),'.tif'],'tif');
%         imwrite(img(end/2+1:end,end/2+1:end),[path_target,'\fullImage_4_',num2str(iE),'.tif'],'tif');
    end
    imgR = img(1:reduction:end,1:reduction:end);
    imwrite(imgR,[path_target,'\img_reduce_res_iE=',num2str(iE),'.tif'],'tif');
end
clock


