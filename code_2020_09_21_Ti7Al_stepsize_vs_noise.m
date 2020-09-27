% The reviewer comment that the step size was too small.
% I would like to confirm my impression, which is that smaller step size
% can result in higher [ratio] of correlated points.
addChenFunction;
% data with step size = 2
dir_2 = 'E:\Ti7Al_E1_insitu_tension\Analysis_selected_area\r6c5';
% data with step size = 5
dir_5 = 'E:\Ti7Al_E1_insitu_tension\SEM Data\RC\r6c5';
% file names are the same
f = 'Ti7Al_E1_S7_r6c5.mat';

%%
d2 = load(fullfile(dir_2,f))
d5 = load(fullfile(dir_5,f))

%% compare percentage of data points with sigma == -1

sum(d2.sigma(:) == -1)/numel(d2.sigma(:))
sum(d5.sigma(:) == -1)/numel(d5.sigma(:))

%% average sigma value. sigma=0 indicates perfect match
sigma_2 = d2.sigma;
sigma_2 = sigma_2(100:end-100,100:end-100);
sigma_2(sigma_2 == -1) = nan;
nanmean(sigma_2(:))

sigma_5 = d5.sigma;
sigma_5 = sigma_5(40:end-40,40:end-40);
sigma_5(sigma_5 == -1) = nan;
nanmean(sigma_5(:))
