% std of DIC data for a pair of calibration images

clear;
clc;

working_dir = 'E:\Ti7Al_E1_insitu_tension\SEM Data\calibration';
% f_ref = 'Ti7Al_Ref_r0c0.mat';
f_def = 'Ti7Al_Ref_r0c5.mat';

load(fullfile(working_dir,f_def));

nanstd(exx(:))  % = 0.0069
nanstd(u(:))    % = 0.2490
nanstd(v(:))    % = 0.1769
