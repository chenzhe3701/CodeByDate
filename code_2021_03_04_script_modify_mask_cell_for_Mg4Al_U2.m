% a script to modify mask_cell
iEs = 0:13;
iBs = 1:14;
inds_to_remove = [1,4,7,4,4, 7,4,7,5,6, 6,6,6,5]; % for each iE, the index of the grain ID to be removed from the analysis of adding gb. 

%%
load('Mg4Al_U2_mask_cell.mat')
old_mask_cell = mask_cell;
mask_cell = cell(1,14);

for iE = 0:13
    iB = iE+1;
    ic = 1;
    for ii=1:length(old_mask_cell{iB})
        if ii~=inds_to_remove(iB)
            mask_cell{iB}{ic} = old_mask_cell{iB}{ii};
            ic = ic+1;
        end
    end
end

%%
save('Mg4Al_U2_mask_cell.mat','mask_cell')