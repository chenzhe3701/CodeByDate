
%% Using Mg4Al_U1 (non-symmetric sample) data, study how to process EBSD data by MTEX

%% install mtex
addChenFunction;
cwd = pwd();
cd('D:\p\m\mtex');
startup_mtex();
cd(cwd);

%% load and plot IPF to illustrate
cd('E:\zhec umich Drive\2020-10-13 Mg4Al insitu EBSD')
file_name = 'Mg4Al 2020-10-13 iE=2.osc';

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

% crystal symmetry. Note this is different from crystalSymmetry.load('Mg-Magnesium.cif') 
CS = crystalSymmetry('6/mmm', [3.2 3.2 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Mg');

% can import and export as ang
ebsd = EBSD.load(file_name, CS);
ebsd.export_ang([replace(file_name,'.osc',''),' export.ang']);

% UMich data, so setting 1.
ebsd = EBSD.load(file_name, CS, 'convertEuler2SpatialReferenceFrame','setting 1');

figure;
ebsd.plot(ebsd.orientations)

%% Explicitly describe the ipfKey. --> If set ipf direction to vector3d.X, we can see the difference between setting 1 and 3  
ipfKey = ipfTSLKey(CS);
ipfKey.inversePoleFigureDirection = vector3d.Z;    % which direction to plot for IPF
ipf_colors = ipfKey.orientation2color(ebsd('Mg').orientations);
figure;
plot(ebsd, ipf_colors);

%% grains
[grains, ebsd.grainId] = calcGrains(ebsd('Mg'),'threshold',10*degree);

figure;
plot(grains.boundary,'linewidth',2)



%% (discrete) pole figures and inverse pole figures
% the fixed crystal direction
h = Miller({0,0,0,1},{1,0,-1,0},ebsd('Mg').CS);
figure;
plotPDF(ebsd.orientations, h, 'figsize','medium');
figure;
plotIPDF(ebsd.orientations, vector3d.Z, 'markersize', 5);

%% crystal shape
cS = crystalShape.hex(ebsd('Mg').CS)
figure;
plot(cS);

%% uninstall mtex
cd('D:\p\m\mtex');
uninstall_mtex();
cd(cwd);