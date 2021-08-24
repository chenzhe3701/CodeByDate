% file_name = cell for multiple files

function [] = pole_figure_by_mtex_multi_file(file_name_cell, varargin)

p = inputParser;
addRequired(p,'file_name_cell');   
addParameter(p,'setting_str','setting 1');
addParameter(p,'clims',[]);
addParameter(p,'axis_label',{'ED','RD'});
parse(p, file_name_cell, varargin{:});

file_name_cell = p.Results.file_name_cell;
setting_str = p.Results.setting_str;
clims = p.Results.clims;
axis_label = p.Results.axis_label;

if isempty(setting_str)
    setting_str = 'setting 1';
    warning('no input, use setting 1');
end

if ~exist('setting_str','var')
    setting_str = 'setting 1';
    warning('no input, use setting 1');
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
       a(ii).String = axis_label{1};
   elseif strcmpi(a(ii).String,'Y')
       a(ii).String = axis_label{2};
   end
end

% change colorbar    
if ~isempty(clims)
    a = findobj(gcf,'Type','axes');
    
    % sort according to location, from left to right
    for ii = 1:length(a)
        axis_pos(ii) = a(ii).Position(1);
    end
    [~, ind] = sort(axis_pos);
    
    for ii=1:length(a)
        axes(a(ind(ii)));  % make current axis as a(ind(ii))
        caxis(clims{ii});
    end
end

%% uninstall mtex
% cd('D:\p\m\mtex-5.3.1');
% uninstall_mtex();
% cd(temp_dir);

end