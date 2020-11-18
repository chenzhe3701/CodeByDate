% analyze the EBSD, pole figures, for 3 areas, each 1200x1200um, for Mg4Al_C2 sample

%% install mtex. 5.4.0 has some problem with calcDensity(ori)
addChenFunction;
temp_dir = pwd();
cd('D:\p\m\mtex-5.3.1');
startup_mtex();
cd(temp_dir);

%%
working_dir = 'E:\zhec umich Drive\2020-11-12 Mg4Al_C2 EBSD big AOI\';
output_dir = 'E:\zhec umich Drive\2020-11-12 Mg4Al_C2 EBSD big AOI\analysis\';

% osc_name = 'Mg4Al_C2_left.osc';
% osc_name = 'Mg4Al_C2_mid.osc';
% osc_name = 'Mg4Al_C2_right.osc';

pos_str = {'left','mid','right'};
for iP = 1:3
    osc_name = ['Mg4Al_C2_', pos_str{iP}, '.osc'];
    
    setMTEXpref('xAxisDirection','east');
    setMTEXpref('zAxisDirection','intoPlane');
    
    CS = crystalSymmetry('6/mmm', [3.2, 3.2, 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Mg');
    ss = specimenSymmetry();
    
    ebsd = EBSD.load(osc_name, CS, 'convertEuler2SpatialReferenceFrame', 'setting 1');
    
    ori = ebsd.orientations;
    % plot(ebsd,ori); % IPF map
    
    odf = calcDensity(ori);
    rot = rotation.byAxisAngle(vector3d.Y, 90*degree);  % rotation to make x out of plan, z to east
    
    % (1) Pole figures of sample normal direction
    figure;
    plotPDF(odf, Miller({0,0,0,1},{1,0,-1,0},{1,1,-2,0},CS), 'antipodal');
    
    mtexColorbar;
    pos = get(gcf,'position');
    pos(3) = 1200; % new width
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
    
    % (2) Pole figures of sample loading direction
    figure;
    plotPDF(rot * odf, Miller({0,0,0,1},{1,0,-1,0},{1,1,-2,0},CS), 'antipodal');
    
    mtexColorbar;
    pos = get(gcf,'position');
    pos(3) = 1200; % new width
    pos(4) = 375; % new height
    set(gcf,'position',pos);
    
    % change label
    a = findall(gcf,'Type','text');
    for ii = 1:length(a)
        if strcmpi(a(ii).String,'X')
            a(ii).String = 'RD';
        elseif strcmpi(a(ii).String,'Y')
            a(ii).String = 'RD';
        end
    end
end
%% Part 2, read grain files, rotate, save as matrix
pos_str = {'left','mid','right'};

for iP = 1:3
    file_name_1 = ['Mg4Al_C2_', pos_str{iP}, ' raw grain_file_type_1.txt'];
    file_name_2 = ['Mg4Al_C2_', pos_str{iP}, ' raw grain_file_type_2.txt'];
    [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,file_name_1));
    [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,file_name_2));
    
    % find column
    column_index_1 = find_variable_column_from_grain_file_header(EBSD_header_1,...
        {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
    column_index_2 = find_variable_column_from_grain_file_header(EBSD_header_2,...
        {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});
    
    % read type-2 grain file and get average info for grains
    gID = EBSD_data_2(:,column_index_2(1));
    gPhi1 = EBSD_data_2(:,column_index_2(2));
    gPhi = EBSD_data_2(:,column_index_2(3));
    gPhi2 = EBSD_data_2(:,column_index_2(4));
    
    % EBSD data, from type-1 grain file. (column, data) pair:
    % (1,phi1) (2,phi) (3,phi2) (4,xMicron) (5,yMicron) (6,IQ) (7,CI) (8,Fit) (9,grain-Id) (10,edgeGrain?)
    % Read EBSD data.  IQ,CI,Fit are not needed for now, but might need in future
    x_um = EBSD_data_1(:,column_index_1(5));
    y_um = EBSD_data_1(:,column_index_1(6));
    unique_x = unique(x_um(:));
    ebsdStepSize = unique_x(2) - unique_x(1);
    mResize = (max(x_um(:)) - min(x_um(:)))/ebsdStepSize + 1;
    nResize = (max(y_um(:)) - min(y_um(:)))/ebsdStepSize + 1;
    
    phi1 = reshape(EBSD_data_1(:,column_index_1(2)),mResize,nResize)';
    phi = reshape(EBSD_data_1(:,column_index_1(3)),mResize,nResize)';
    phi2 = reshape(EBSD_data_1(:,column_index_1(4)),mResize,nResize)';
    % change it to degrees, if necessary
    if max(phi1(:))<7 && max(phi(:))<7 && max(phi2(:))<7
        phi1 = phi1*180/pi();
        phi = phi*180/pi();
        phi2 = phi2* 180/pi();
    end
    x_um = reshape(EBSD_data_1(:,column_index_1(5)),mResize,nResize)';
    y_um = reshape(EBSD_data_1(:,column_index_1(6)),mResize,nResize)';
    ID = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
    
    boundary = find_boundary_from_ID_matrix(ID);
    
    [phi1, phi, phi2] = align_euler_to_sample(phi1, phi, phi2, 'custom', 90, 180, 0);
    [gPhi1, gPhi, gPhi2] = align_euler_to_sample(gPhi1, gPhi, gPhi2, 'custom', 90, 180, 0);
    euler_rotated = true;
    
    save(fullfile(working_dir, ['Mg4Al_C2_', pos_str{iP}, '_EBSD.mat']), 'x_um','y_um', 'ID','phi1','phi','phi2','gID','gPhi1','gPhi','gPhi2','euler_rotated');
end
%% check IPF color map
% [nR,nC] = size(phi1);
% IPF_map = zeros(nR,nC,3);
% for iR = 1:nR
%     for iC = 1:nC
%         color = calculate_IPF_color_hcp([phi1(iR,iC), phi(iR,iC), phi2(iR,iC)], [0,0,0], [0 0 1]);
%         IPF_map(iR,iC,:) = reshape(color,1,1,3);
%     end
% end
%% uninstall mtex
cd('D:\p\m\mtex-5.3.1');
uninstall_mtex();
