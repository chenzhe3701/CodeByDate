% 2020-12-30 make movie using .tif images, with the selected prefix, from the selected folder 
% required parameters:
% 'folder'
% 'img_prefix'
% 'img_suffix'

function [] = make_movie_using_images(varargin)

p = inputParser;

p.addParameter('folder',[]);
p.addParameter('img_prefix',[]);
p.addParameter('img_suffix',[]);

p.parse(varargin{:});

folder = p.Results.folder;
img_prefix = p.Results.img_prefix;
img_suffix = p.Results.img_suffix;

outputVideo = VideoWriter(fullfile(folder,'output_video.avi'));
outputVideo.FrameRate = 1;  % frame per second
open(outputVideo)

fileList = dir([folder,'\',img_prefix,'*',img_suffix]);

folder_name = fileList(1).folder;
file_names = sort({fileList.name});
[file_names,index] = sort_nat(file_names);

for ii = 1:length(fileList)
    i_img = imread(fullfile(folder_name,file_names{ii}));
    i_frame(ii) = im2frame(i_img);
    writeVideo(outputVideo,i_img);    
end
close(outputVideo);

end