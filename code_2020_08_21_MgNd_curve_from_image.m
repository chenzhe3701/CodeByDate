% get stress strain from image, for Mg2.4Nd
% 2020-08-21

%% open image, crop curve area
I = imread('data\zhihua_MgNd_tensile_stressstrain.png');
I = imcrop(I, [82, 48, 440-82, 242-48]);

figure; 
image(I);
%% draw a line on top of the curve
H = drawfreehand;

%% make a mask
mask = H.createMask;

%% save the mask
save('data\MgNd curve mask.mat','mask');

%%  process
% find curve from mask
mask = load('data\MgNd curve mask.mat','mask');
mask = mask.mask;
mask = imresize(mask, [181,141]);

[nr,nc] = size(mask);
curve = zeros(nr,nc);
for ic = 1:nc
    ir = find(mask(:,ic), 1, 'first');
    curve(ir,ic) = 1;
end

% convert to stress and strain
curve = flipud(curve);
[nr,nc] = size(curve);
stress = 0;
strain = 0;
iL =1;
for ic = 1:nc
    ir = find(curve(:,ic), 1, 'first');
    stress = [stress; ir];
    strain = [strain; ic/1000];
end

% plot curve
figure;
plot(strain, stress);

%% save
save('data\Mg2.4Nd tensile from image.mat','stress','strain');