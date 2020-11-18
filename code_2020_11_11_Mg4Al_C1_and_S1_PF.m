% Analyze the texture of Mg4Al samples
% chenzhe, 2020-11-11

%% install mtex. 5.4.0 has some problem with calcDensity(ori)
addChenFunction;
cd('D:\p\m\mtex-5.3.1');
startup_mtex();

%% saving dir
saving_dir = 'C:\Users\ZheChen\Desktop\Mg4Al Comparison';

%% Mg4Al_C1 --> note setting-1
% working dir for Mg4Al_C1
working_dir_Mg4Al_C1 = 'E:\zhec umich Drive\2020-10-23 Mg4Al_C1 insitu EBSD';
file_name = 'Mg4Al_C1_iE=0.osc';
cd(working_dir_Mg4Al_C1);

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

% crystal symmetry. Note this is different from crystalSymmetry.load('Mg-Magnesium.cif') 
CS = crystalSymmetry('6/mmm', [3.2 3.2 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Mg');     

% UMich data, so setting 1.
ebsd = EBSD.load(file_name, CS, 'convertEuler2SpatialReferenceFrame','setting 1');

ori = ebsd.orientations;
plot(ebsd,ori);
print(fullfile(saving_dir,'Mg4Al_C1_IPF.tiff'),'-dtiff');

odf = calcDensity(ori);
rot = rotation.byAxisAngle(vector3d.Y, 90*degree);  % rotation to make x out of plan, z to east
figure;
plotPDF(odf, Miller({0,0,0,1},{1,0,-1,0},CS), 'antipodal');
mtexColorbar;
pos = get(gcf,'position');
pos(3) = 800; % new width
pos(4) = 375; % new height
set(gcf,'position',pos);
% change label
a = findall(gcf,'Type','text');
for ii = 1:length(a)
   if strcmpi(a(ii).String,'X')
       a(ii).String = 'ED';
   elseif strcmpi(a(ii).String,'Y')
       a(ii).String = 'RD';
   end
end

print(fullfile(saving_dir,'Mg4Al_C1_PF.tiff'),'-dtiff');
%% Mg4Al_S1 --> note setting-2
% working dir for Mg4Al_C1
working_dir_Mg4Al_S1 = 'E:\Mg4Al_S1_insitu\EBSD Data';
file_name = 'Mg4Al s1.osc';
cd(working_dir_Mg4Al_S1);

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

% crystal symmetry. Note this is different from crystalSymmetry.load('Mg-Magnesium.cif') 
CS = crystalSymmetry('6/mmm', [3.2 3.2 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Mg');     

% UCSB data, so setting 1.
ebsd = EBSD.load(file_name, CS, 'convertEuler2SpatialReferenceFrame','setting 2');
ebsd = ebsd(ebsd.inpolygon([100,80,600,660]));  % --> select the approximate area 

ori = ebsd.orientations;
plot(ebsd,ori);
print(fullfile(saving_dir,'Mg4Al_S1_IPF.tiff'),'-dtiff');

odf = calcDensity(ori);
figure;
plotPDF(odf, Miller({0,0,0,1},{1,0,-1,0},CS), 'antipodal');
mtexColorbar;
pos = get(gcf,'position');
pos(3) = 800; % new width
pos(4) = 375; % new height
set(gcf,'position',pos);
% change label
a = findall(gcf,'Type','text');
for ii = 1:length(a)
   if strcmpi(a(ii).String,'X')
       a(ii).String = 'ED';
   elseif strcmpi(a(ii).String,'Y')
       a(ii).String = 'RD';
   end
end

print(fullfile(saving_dir,'Mg4Al_S1_PF.tiff'),'-dtiff');
%% uninstall mtex
cd('D:\p\m\mtex-5.3.1');
uninstall_mtex();
