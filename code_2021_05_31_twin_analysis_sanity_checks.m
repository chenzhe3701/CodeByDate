%% for the Mg4Al_U1 persistent twin, twin-retwin-detwin analysis

%% [1, good] first, each past_or_present_twin_grain_ID should only have 1 label.
for iE = 0:13
    iB = iE + 1;
    IDs = unique(past_or_present_twin_grain_ID_cell{iB});
    for ii = 1:length(IDs)
        ID_current = IDs(ii);
        inds = ismember(past_or_present_twin_grain_ID_cell{iB}, ID_current);
        labels = unique(past_or_present_twin_grain_label_cell{iB}(inds));
        if length(labels)>1
            error('aa');
        end
    end
end

disp('OK 1');
    
%% [2] Each 'current twin grain' should only belong to one 'past_or_present_twin grain', except also overlap with '0'

for iE = 0:13
    iB = iE + 1;
    gList_current = nan_unique(ID_c_iB_to_1_cell{iB});
    gList_current(isnan(gList_current)) = [];
    
    for ii = 1:length(gList_current)
        ID_current = gList_current(ii);
        
        inds = ismember(ID_c_iB_to_1_cell{iB}, ID_current);
        % determine if twinned
        if any(variant_pixel_iB_to_1_cell{iB}(inds))
            ids_past_or_present = unique(past_or_present_twin_grain_ID_cell{iB}(inds));
            ids_past_or_present(ids_past_or_present==0) = [];
            if length(ids_past_or_present) > 1
                error('aa');
            end
        end
    end
end

disp('OK 2');

%% [3] the non overlap area should not have any past_or_present_twin_grain info
