% process the external image taken on 2020-12-01
% image 'r3 s10.tiff' was scanned vertically, but stored horizontally.
% We can flip it along the diagonal:upperLeft-lowerRight to make it show
% the real object.

working_dir = 'E:\zhec umich Drive\2020-12-01 image with external scan';

I = imread('test_image r3 s10.tiff');

% For this image, the scan direction was vertical: for each col left->right
% (x: 0->max), scan each row top->bottom (y: 0->max).
% It was saved as each scan line as a horizontal line.
% Therefore, we can rotate 90 degree ccw, then flipud, so that the image shows the real object.

I = imrotate(I, 90);
I = flipud(I);
imwrite(I, 'test_image r3 s10 realOjb verticalScan.tiff');