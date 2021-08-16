% file_name = cell for multiple files

function [] = pole_figure_by_mtex_multi_file(file_name_cell, setting_str, max_intensity)

if ~exist('setting_str','var')
    setting_str = 'setting 1';
    warning('no input, use setting 1');
end

if ~exist('file_name_cell','var')
    [file_name_cell,path_name] = uigetfile('E:\zhec umich Drive');
    file_name_cell = fullfile(path_name,file_name_cell);
end

temp_dir = pwd();
cd('D:\p\m\mtex-5.3.1');
startup_mtex();
cd(temp_dir);

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

% crystal symmetry. Note this is different from crystalSymmetry.load('Mg-Magnesium.cif') 
CS = crystalSymmetry('6/mmm', [3.2 3.2 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Mg');     

% can use a cell to include multiple files
ebsd = [];
for ii = 1:length(file_name_cell)
    % UMich data, so setting 1.
    t = EBSD.load(file_name_cell{ii}, CS, 'convertEuler2SpatialReferenceFrame', setting_str);
    if ii == 1
        ebsd = t;
    else
        dd = t(2).x -t(1).x;
        t.x = t.x + max(ebsd.x) + dd;
        ebsd = [ebsd,t];
    end
end

    
ori = ebsd.orientations;

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

if exist('max_intensity','var')
   a = findobj(gcf,'Type','axes');
    for ii = 1:length(a)
        axis_pos(ii) = a(ii).Position(1);
    end
    [~, ind] = sort(axis_pos);
    
    for ii=1:length(a)
         axes(a(ind(ii)));  % make current axis as a(ind(ii))
         caxis([0 max_intensity(ii)]);
    end
end

%% uninstall mtex
% cd('D:\p\m\mtex-5.3.1');
% uninstall_mtex();
% cd(temp_dir);

end