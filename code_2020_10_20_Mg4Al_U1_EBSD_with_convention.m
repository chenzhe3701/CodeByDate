% Performed EBSD on a tested Mg4Al sample which was stored in air for 10 days. 
% This script include some discussion on the MTEX vs OIM convention


%% install mtex
addChenFunction;
data_path = 'E:\zhec umich Drive\2020-10-19 EBSD quality';

cd('D:\p\m\mtex');
startup_mtex();
cd(data_path);

%% load data, generate IPF, export .ang
plotx2east;
plotzIntoPlane;
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

osc_file_name = '20kV BI=16 4x4 binning FPS=33.osc'; 

CS_Mg = crystalSymmetry('6/mmm', [3.2, 3.2, 5.2], 'x||a', 'mineral','Mg');
CS_Mg17Al12 = crystalSymmetry('-43m', [11 11 11], 'mineral','Mg17Al12');
CS = {'notIndexed', CS_Mg, CS_Mg17Al12};

ebsd = EBSD.load(osc_file_name, CS, 'convertEuler2SpatialReferenceFrame', 'setting 1');
ebsd.export_ang('testexport 4x4.ang');

%% MTEX's crystal symmetry file seems to have different convention from that used by TSL OIM analysis
osc_file_name = '20kV BI=16 1x1 binning FPS=33.osc'; 

% (1) Simply load .osc, and output .ang, the CS is '622', which CANNOT be understood by OIM analysis 
ebsd = EBSD.load(osc_file_name);
ebsd.export_ang('[1] 1x1 FPS=33 osc_to_ang 622.ang');

% (2) If load Mg crystal symmetry from the .cif file that came with MTEX, there is a 30 deg rotation 
CS = {'notIndexed', crystalSymmetry('6/mmm', [3.2 3.2 5.2], 'X||a*', 'Y||b', 'Z||c', 'mineral', 'Mg')};
ebsd = EBSD.load(osc_file_name, CS, 'interface','osc');
ebsd.export_ang('[2] 1x1 FPS=33 osc_to_ang 6mmm HKL.ang');

% (3) Correct for the 30 deg rotation 
CS = {'notIndexed', crystalSymmetry('6/mmm', [3.2 3.2 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Mg')};
ebsd = EBSD.load(osc_file_name, CS, 'interface','osc');
ebsd.export_ang('[3] 1x1 FPS=33 osc_to_ang 6mmm TSL.ang');

% (4) If load Mg crystal symmetry from the .cif file that came with MTEX, there is a 30 deg rotation 
CS = {'notIndexed', crystalSymmetry('62', [3.2 3.2 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Mg')};
ebsd = EBSD.load(osc_file_name, CS, 'interface','osc');
ebsd.export_ang('[4] 1x1 FPS=33 osc_to_ang 62 TSL.ang');
%% Note: Even defined as '62', it still export as '622', so need to manually correct to '62', or modify the function ...


%% Additionally, the crystal symmetry definition TSL 'X||a' vs HKL 'X||a*' affects crystal lattice illustration of grain orientation. 
% Therefore, it's better to correct to TSL 'X||a'.

% [5] HKL definition gives a different crystal lattice illustration
CS = {'notIndexed', crystalSymmetry('6/mmm', [3.2 3.2 5.2], 'X||a*', 'Y||b', 'Z||c', 'mineral', 'Mg')};
ebsd = EBSD.load(osc_file_name, CS, 'interface','osc','convertEuler2SpatialReferenceFrame');
disp('HKL definition, convertEuler2SpatialReferenceFrame:');
ebsd.orientations(1)
cS = crystalShape.hex(ebsd.CS)
figure;
scaling = 30;
ebsd.plot(ebsd.prop.confidenceindex); 
hold on;
plot(64,64,-50, ebsd(64,64).orientations * cS * scaling);

% [6] Correct definition
CS = {'notIndexed', crystalSymmetry('6/mmm', [3.2 3.2 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Mg')};
ebsd = EBSD.load(osc_file_name, CS, 'interface','osc','convertEuler2SpatialReferenceFrame');
disp('TSL definition, convertEuler2SpatialReferenceFrame:');
ebsd.orientations(1)
cS = crystalShape.hex(CS{2})
figure; 
ebsd.plot(ebsd.prop.confidenceindex);
hold on;
plot(64,64,-50, ebsd(64,64).orientations * cS * scaling);

% If just plot crystal, without IPF overlay, need to correct the view angle, as MTEX default to view(3)
figure;
plot(ebsd(64,64).orientations * cS);
view(0,-90);    % --> remember to correct for the view angle. MTEX uses default 3D view, which gives the wrong crystal orientation on IPF map.
hold on;
plot3([0,1],[0,0],[0,0],'-r','linewidth',3)
plot3([0,0],[0,1],[0,0],'-g','linewidth',3)
plot3([0,0],[0,0],[0,1],'-b','linewidth',3)

% If forget to set the correct view angle
figure;
plot(ebsd(64,64).orientations * cS);
hold on;
plot3([0,1],[0,0],[0,0],'-r','linewidth',3)
plot3([0,0],[0,1],[0,0],'-g','linewidth',3)
plot3([0,0],[0,0],[0,1],'-b','linewidth',3)



%%
ori = orientation.byEuler(0*degree,0*degree,0*degree,CS{2})
cS = crystalShape.hex(ebsd.CS)
plot(cS)
xlabel('x')
zlabel('z')
axis on;

%% uninstall mtex
cd('D:\p\m\mtex');
uninstall_mtex();
cd(cwd);