%% Analyze the twin thickness
% Also, it's better to create a function for segmenting twin
addChenFunction;
%% [1] Look at Mg4Al_U2, resolution = 1um/dp
working_dir = 'E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD';
sample_name = 'Mg4Al_U2';
output_dir = fullfile(working_dir, 'analysis', 'twin thickness analysis');
mkdir(output_dir);

strain_sg = [0, -0.0075, -0.015, -0.025, ...
    -0.023, -0.017, -0.0075, 0.002, ...
    -0.0075, -0.015, -0.025, ...
    -0.017, -0.0075, 0];

% strain values, Mg4Al_U2
strain_ebsd = [0, -0.0085, -0.0220, -0.0342, ...
    -0.0343, -0.0268, -0.0150, -0.0032, ...
    -0.0112, -0.0207, -0.0297, ...
    -0.0235, -0.0130, -0.0045];

% load twin variant map. ==> this should be corrected later. iB vs iE issue.  
d = load(fullfile(working_dir, 'analysis', 'variant_maps.mat'));
variant_point_wise = d.variant_point_wise;
if length(variant_point_wise)==13
    variant_point_wise(2:14) = variant_point_wise(1:13);
    variant_point_wise{1} = zeros(size(variant_point_wise{1}));
end

% Load data directly from regular grain file. We currently can only analyze
% grains in the unmodified grain file, because the fitted ellipse
% information was not processed previously.
% For each child grain, check the variant map to see if it is a twin.

for iE = 1:13
    iB = iE + 1;
    
    data = grain_file_to_data(fullfile(working_dir,[sample_name,' grain_file_type_1 iE=',num2str(iE),'.txt']), ...
        fullfile(working_dir,[sample_name,' grain_file_type_2 iE=',num2str(iE),'.txt']));
    
    ID = data.ID;
    phi1 = data.phi1;
    phi = data.phi;
    phi2 = data.phi2;
    x = data.x;
    y = data.y;
    
    gID = data.gID;
    gCenterX = data.gCenterX;
    gCenterY = data.gCenterY;
    gAspectRatio = data.gAspectRatio;
    gMajorAxis = data.gMajorAxis;
    gMinorAxis = data.gMinorAxis;
    gMajorOrientation = data.gMajorOrientation;
    gMaxFeretDia = data.gMaxFeretDia;
    gMinFeretDia = data.gMinFeretDia;
    
    variant_map = variant_point_wise{iB};
    boundary = find_one_boundary_from_ID_matrix(ID);
    
    myplot(x,y,variant_map,boundary);
    make_variant_map_background('bg',[1,1,1],'tickOff',false);
        
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    title(['load step ',num2str(iE),', \epsilon=',num2str(strain_sg(iB),'%.3f')], 'fontweight','normal');
    set(gca,'fontsize',14);
    
    thickness = [];
    for ii = 1:length(gID)
        ID_current = gID(ii);
        
        inds = ismember(ID, ID_current);
        isTwin = mode(variant_map(inds))>0; % mode of child grain's twin variant label
        
        if isTwin
            ind = find(gID == ID_current);
            gx = gCenterX(ind);
            gy = gCenterY(ind);
            
            a = gMajorAxis(ind);
            b = gMinorAxis(ind);
            alpha = gMajorOrientation(ind);
            
            thickness = [thickness; b];
            
            plot3_ellipse(a,b,gx,gy,alpha);
        end
        
        thickness_cell{iB} = thickness * 2;
    end
    
    print(fullfile(output_dir, ['twin thick iE=',num2str(iE),'.tiff']), '-dtiff');
    
    close all;
end

save(fullfile(output_dir, 'twin_thick_data.mat'), 'thickness_cell');

close all;
% data, and grouping variable
t = [];
g = [];
for iE = 1:13
    iB = iE + 1;
    t = [t; thickness_cell{iB}];
    g = [g; repmat(iE,length(thickness_cell{iB}),1)];
end
figure;
boxplot(t,g);
xlabel('Load step');
ylabel('Thickness (\mum)');
set(gca,'fontsize',14);
print(fullfile(output_dir, 'twin thick box plot.tiff'), '-dtiff');

%% [2] Look at UM134_Mg_C3, resolution = 1um/dp
working_dir = 'E:\zhec umich Drive\2021-01-29 UM134 Mg_C3 insitu EBSD';
sample_name = 'UM134_Mg_C3';
output_dir = fullfile(working_dir, 'analysis', 'twin thickness analysis');
mkdir(output_dir);

strain_sg = [0, -0.0075, -0.015, -0.025, ...
    -0.023, -0.017, -0.0075, 0.002, ...
    -0.0075, -0.015, -0.025, ...
    -0.017, -0.0075, 0];

% strain values, UM134_Mg_C3
strain_ebsd = [0.0000, -0.0292, -0.0385, -0.0465, ... 
    -0.0414, -0.0313, -0.0145, -0.0012, ...
    -0.0179, -0.0320, -0.0438, ...
    -0.0301, -0.0145, -0.0016];

% load twin variant map. ==> this should be corrected later. iB vs iE issue.  
d = load(fullfile(working_dir, 'analysis', 'variant_maps.mat'));
variant_point_wise = d.variant_point_wise;
if length(variant_point_wise)==13
    variant_point_wise(2:14) = variant_point_wise(1:13);
    variant_point_wise{1} = zeros(size(variant_point_wise{1}));
end

% Load data directly from regular grain file. We currently can only analyze
% grains in the unmodified grain file, because the fitted ellipse
% information was not processed previously.
% For each child grain, check the variant map to see if it is a twin.

for iE = 1:13
    iB = iE + 1;
    
    data = grain_file_to_data(fullfile(working_dir,[sample_name,' grain_file_type_1 iE=',num2str(iE),'.txt']), ...
        fullfile(working_dir,[sample_name,' grain_file_type_2 iE=',num2str(iE),'.txt']));
    
    ID = data.ID;
    phi1 = data.phi1;
    phi = data.phi;
    phi2 = data.phi2;
    x = data.x;
    y = data.y;
    
    gID = data.gID;
    gCenterX = data.gCenterX;
    gCenterY = data.gCenterY;
    gAspectRatio = data.gAspectRatio;
    gMajorAxis = data.gMajorAxis;
    gMinorAxis = data.gMinorAxis;
    gMajorOrientation = data.gMajorOrientation;
    gMaxFeretDia = data.gMaxFeretDia;
    gMinFeretDia = data.gMinFeretDia;
    
    variant_map = variant_point_wise{iB};
    boundary = find_one_boundary_from_ID_matrix(ID);
    
    myplot(x,y,variant_map,boundary);
    make_variant_map_background('bg',[1,1,1],'tickOff',false);
        
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    title(['load step ',num2str(iE),', \epsilon=',num2str(strain_sg(iB),'%.3f')], 'fontweight','normal');
    set(gca,'fontsize',14);
    
    thickness = [];
    for ii = 1:length(gID)
        ID_current = gID(ii);
        
        inds = ismember(ID, ID_current);
        isTwin = mode(variant_map(inds))>0; % mode of child grain's twin variant label
        
        if isTwin
            ind = find(gID == ID_current);
            gx = gCenterX(ind);
            gy = gCenterY(ind);
            
            a = gMajorAxis(ind);
            b = gMinorAxis(ind);
            alpha = gMajorOrientation(ind);
            
            thickness = [thickness; b];
            
            plot3_ellipse(a,b,gx,gy,alpha);
        end
        
        thickness_cell{iB} = thickness * 2;
    end
    
    print(fullfile(output_dir, ['twin thick iE=',num2str(iE),'.tiff']), '-dtiff');
    
    close all;
end

save(fullfile(output_dir, 'twin_thick_data.mat'), 'thickness_cell');

close all;
% data, and grouping variable
t = [];
g = [];
for iE = 1:13
    iB = iE + 1;
    t = [t; thickness_cell{iB}];
    g = [g; repmat(iE,length(thickness_cell{iB}),1)];
end
figure;
boxplot(t,g);
xlabel('Load step');
ylabel('Thickness (\mum)');
set(gca,'fontsize',14);
print(fullfile(output_dir, 'twin thick box plot.tiff'), '-dtiff');



%% [3] UM129_Mg_C2 data, resolution = 5 um/dp
working_dir = 'E:\zhec umich Drive\2021-08-20 UM129 Mg_C2 insitu EBSD';
sample_name = 'UM129_Mg_C2';

output_dir = fullfile(working_dir, 'analysis', 'twin thickness analysis');
mkdir(output_dir);

% strain value, UM129_Mg_C2
strain_ebsd = [0.0000, -0.0110, -0.0179, -0.0277, ...
    -0.0254, -0.0201, -0.0102, -0.0024, ...
    -0.0100, -0.0181, -0.0281, ...
    -0.0204, -0.0104, -0.0030];

% load twin variant map. ==> this should be corrected later. iB vs iE issue.  
d = load(fullfile(working_dir, 'analysis', 'variant_maps.mat'));
variant_point_wise = d.variant_point_wise;

% Load data directly from regular grain file. We currently can only analyze
% grains in the unmodified grain file, because the fitted ellipse
% information was not processed previously.
% For each child grain, check the variant map to see if it is a twin.

for iE = 0:13
    iB = iE + 1;
    
    data = grain_file_to_data(fullfile(working_dir,[sample_name,' grain_file_type_1 iE=',num2str(iE),'.txt']), ...
        fullfile(working_dir,[sample_name,' grain_file_type_2 iE=',num2str(iE),'.txt']));
    
    ID = data.ID;
    phi1 = data.phi1;
    phi = data.phi;
    phi2 = data.phi2;
    x = data.x;
    y = data.y;
    
    gID = data.gID;
    gCenterX = data.gCenterX;
    gCenterY = data.gCenterY;
    gAspectRatio = data.gAspectRatio;
    gMajorAxis = data.gMajorAxis;
    gMinorAxis = data.gMinorAxis;
    gMajorOrientation = data.gMajorOrientation;
    gMaxFeretDia = data.gMaxFeretDia;
    gMinFeretDia = data.gMinFeretDia;
    
    variant_map = variant_point_wise{iB};
    boundary = find_one_boundary_from_ID_matrix(ID);
    
    myplot(x,y,variant_map,boundary);
    make_variant_map_background('bg',[1,1,1],'tickOff',false);
        
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    title(['load step ',num2str(iE),', \epsilon=',num2str(strain_sg(iB),'%.3f')], 'fontweight','normal');
    set(gca,'fontsize',14);
    
    thickness = [];
    for ii = 1:length(gID)
        ID_current = gID(ii);
        
        inds = ismember(ID, ID_current);
        isTwin = mode(variant_map(inds))>0; % mode of child grain's twin variant label
        
        if isTwin
            ind = find(gID == ID_current);
            gx = gCenterX(ind);
            gy = gCenterY(ind);
            
            a = gMajorAxis(ind);
            b = gMinorAxis(ind);
            alpha = gMajorOrientation(ind);
            
            thickness = [thickness; b];
            
            plot3_ellipse(a,b,gx,gy,alpha);
        end
        
        thickness_cell{iB} = thickness * 2;
    end
    
    print(fullfile(output_dir, ['twin thick iE=',num2str(iE),'.tiff']), '-dtiff');
    
    close all;
end

save(fullfile(output_dir, 'twin_thick_data.mat'), 'thickness_cell');

close all;
% data, and grouping variable
t = [];
g = [];
for iE = 1:13
    iB = iE + 1;
    t = [t; thickness_cell{iB}];
    g = [g; repmat(iE,length(thickness_cell{iB}),1)];
end
figure;
boxplot(t,g);
xlabel('Load step');
ylabel('Thickness (\mum)');
set(gca,'fontsize',14);
print(fullfile(output_dir, 'twin thick box plot.tiff'), '-dtiff');

%% [final] Compare the 3 materials, on 4 plots, each with half cycle
close all;

p = fullfile('E:\zhec umich Drive\2021-02-26 Mg4Al_U2 EBSD','analysis', 'twin thickness analysis');
d = matfile(fullfile(p, 'twin_thick_data.mat'));
c1 = d.thickness_cell;

p = fullfile('E:\zhec umich Drive\2021-01-29 UM134 Mg_C3 insitu EBSD','analysis', 'twin thickness analysis');
d = matfile(fullfile(p, 'twin_thick_data.mat'));
c2 = d.thickness_cell;


p = fullfile('E:\zhec umich Drive\2021-08-20 UM129 Mg_C2 insitu EBSD','analysis', 'twin thickness analysis');
d = matfile(fullfile(p, 'twin_thick_data.mat'));
c3 = d.thickness_cell;

figure('position',[200,200,1400,500]);
for iE = 1:7    
    subplot(1,7,iE);
    iB = iE + 1;
    
    t = [];
    g = [];
    
    t = [t; c1{iB}];
    g = [g; repmat(categorical(cellstr('Mg4Al')), length(c1{iB}), 1)];
    t = [t; c2{iB}];
    g = [g; repmat(categorical(cellstr('Mg (FG)')), length(c2{iB}), 1)];
    
    % t = [t; c3{iB}];
    % g = [g; repmat(categorical(cellstr('UM129_Mg_C2')), length(c3{iB}), 1)];
    
    boxplot(t,g);
    % ylabel('Thickness (\mum)');
    title(['Load step ',num2str(iE)], 'fontweight','normal');
    set(gca,'fontsize',12, 'ylim', [0 55]);

end


figure('position',[200,200,1400,500]);
for iE = 8:13    
    subplot(1,7,iE-7);
    iB = iE + 1;
    
    t = [];
    g = [];
    
    t = [t; c1{iB}];
    g = [g; repmat(categorical(cellstr('Mg4Al')), length(c1{iB}), 1)];
    t = [t; c2{iB}];
    g = [g; repmat(categorical(cellstr('Mg (FG)')), length(c2{iB}), 1)];
    
    % t = [t; c3{iB}];
    % g = [g; repmat(categorical(cellstr('UM129_Mg_C2')), length(c3{iB}), 1)];
    
    boxplot(t,g);
    % ylabel('Thickness (\mum)');
    title(['Load step ',num2str(iE)], 'fontweight','normal');
    set(gca,'fontsize',12, 'ylim', [0 55]);

end












