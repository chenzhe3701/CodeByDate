% codoe_2022_10_02
% organize dirs for wodes working on re-installed laptop

% [data_dir, sample_name, sample_ID, plot_symbol, group_number, ymax]
cells = {'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD\analysis', 'Mg4Al_U2', 'Mg4Al U2', 'o', 1, 300;
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2020-12-23 Mg4Al_C3 insitu EBSD\analysis', 'Mg4Al_C3', 'Mg4Al C3', 's', 1, 300;
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-10-28 Mg4Al_A1 insitu EBSD\analysis', 'Mg4Al_A1', 'Mg4Al A1', 'o', 2, 350;
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-11-05 Mg4Al_A2 insitu EBSD\analysis', 'Mg4Al_A2', 'Mg4Al A2', 's', 2, 350;
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-12-02 Mg4Al_B1 insitu EBSD\analysis', 'Mg4Al_B1', 'Mg4Al B1', 'o', 3, 200;
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-12-04 Mg4Al_B2 insitu EBSD\analysis', 'Mg4Al_B2', 'Mg4Al B2', 's', 3, 200;
    % 'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2020-12-05 UM134_Mg_C1 insitu EBSD\analysis', 'UM134_Mg_C1', 'Mg UM134 C1', 'x', 4, 450;
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-01-15 UM134_Mg_C2 insitu EBSD\analysis', 'UM134_Mg_C2', 'Mg UM134 C2', 'o', 4, 450;
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-01-29 UM134_Mg_C3 insitu EBSD\analysis', 'UM134_Mg_C3', 'Mg UM134 C3', 's', 4, 450;
    % 'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-06-29 UM129_Mg_C1 insitu EBSD\analysis', 'UM129_Mg_C1', 'Mg UM129 C1', 'x', 5, 200;
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-08-20 UM129_Mg_C2 insitu EBSD\analysis', 'UM129_Mg_C2', 'Mg UM129 C2', 'o', 5, 200;
    'H:\Other computers\My Laptop w541 2022\zhec umich Drive\2021-09-03 UM129_Mg_C3 insitu EBSD\analysis', 'UM129_Mg_C3', 'Mg UM129 C3', 's', 5, 200};

%%

for iCell = 1
    winopen(fullfile(cells{iCell,1}, 'twin pct table.tiff'));
end