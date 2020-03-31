% 2020-02-13
% For Ti7Al-E1 sample, collabrated with MSU BYU
% find the position for each FOV, and plot lines to indicate them

addChenFunction;
% load translation data
load('E:\Ti7Al_E1_insitu_tension\Ti7Al_E1_Images\stitched_img\translations_searched_vertical_stop_0.mat');
% load x y pixel coord data
load('E:\Ti7Al_E1_insitu_tension\Analysis_by_Matlab\Ti7Al_E1_EbsdToSemForTraceAnalysis.mat', 'X','Y','boundaryTFB');
% load strain data
load('E:\Ti7Al_E1_insitu_tension\Ti7Al_E1_Images\stitched_DIC\global_method_stitched_7.mat', 'exx');

% plot
myplot(X,Y,exx,boundaryTFB);

for iR = 1:8
   for iC = 1:8
       xt = transX(iR,iC);
       yt = transY(iR,iC);
       imrect(gca,[xt,yt,4096,4096]);
   end
end
