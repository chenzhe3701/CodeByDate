% study the CPFE output for WE43-T6 sample

working_dir = 'E:\Globus download\results WE43-SEM-DIC OneLayerCalibration_SampleSize4_500by400';
file_name = 'QuadratureOutputs16799.csv';
output_dir = 'E:\zhec umich Drive\0_temp_output\WE43-T6 CPFE output analysis';
mkdir(output_dir);

ref_data = 'E:\zhec umich Drive\WE43_T6_C1 EffStrain\EffStrainCombined_iE_1to5.mat';
d_ref = matfile(ref_data);

d = importdata(fullfile(working_dir,file_name),',',0);

%% effective strain reference
eEff = d_ref.EffStrain_fromCPFE_Global;
eEff = eEff{5};

%% x 500, y 400 elements
% col 1:4 = grain ID, x, y, z
% col 5:7 = rot
% col 8:16 = Fe
% col 17:25 = Fp
% col 26:34 = T
% col 35:118 = slip

nR = 400 * 2;
nC = 500 * 2;
nP = 1 * 2;

d = sortrows(d,[4,2,3]);

x = reshape(d(1:end/2,2), nR,nC);
y = reshape(d(1:end/2,3), nR,nC);
ID = reshape(d(1:end/2,1), nR,nC);

%%
for iC = 1:size(d,2)
    myplot(reshape(d(1:end/2, iC), nR,nC));
    print(fullfile(output_dir,['iC=',num2str(iC),'.tiff']), '-dtiff');
    close;
end


