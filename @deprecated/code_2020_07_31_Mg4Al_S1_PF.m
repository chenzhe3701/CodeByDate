%% 2021-06-15 note: we can use the updated code 'D:\p\m\CodeByDate\code_2021_02_07_make_PF_all_samples.m';

%% start mtex
cd('D:\p\m\mtex-5.3.1');
startup_mtex();

%% path to files
pname = 'E:\Mg4Al_S1_insitu\EBSD Data';
% which files to be imported
fname = [pname '\Mg4Al_s1.ang'];
cd(pname);

%%
cs = crystalSymmetry('622', [3.2 3.2 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Magnesium', 'color', [0.53 0.81 0.98]);
CS = {'notIndexed', cs};

setMTEXpref('xAxisDirection','east');      % or 'east', ...
setMTEXpref('zAxisDirection','intoPlane'); % or 'intoPlane'

ebsd = EBSD.load(fname,CS,'interface','ang','convertEuler2SpatialReferenceFrame','setting 2');    
ebsd = ebsd('Magnesium');
ori = ebsd.orientations;
figure; plot(ebsd, ori);

odf = calcDensity(ori);
figure;
plotPDF(odf, Miller({0,0,0,1},{1,0,-1,0},cs), 'antipodal');
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
%% uninstall mtex
cd('D:\p\m\mtex-5.3.1');
uninstall_mtex();

%% print
print('c:\users\ZheChen\Desktop\pf.tiff','-dtiff','-r300') 