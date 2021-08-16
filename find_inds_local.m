% find local indices of grain (ID_target) with ID map
% just replace 5 lines codes with 1 line
% chenzhe, 2021-04-19
function [indR_min, indR_max, indC_min, indC_max] = find_inds_local(ID, ID_target)
    ind_local = ismember(ID, ID_target);
    indC_min = find(sum(ind_local, 1), 1, 'first');
    indC_max = find(sum(ind_local, 1), 1, 'last');
    indR_min = find(sum(ind_local, 2), 1, 'first');
    indR_max = find(sum(ind_local, 2), 1, 'last');
end