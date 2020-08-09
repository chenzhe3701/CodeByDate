% Hough peaks of a sin curve, for study

close all;
addChenFunction;

mat = zeros(200);
for c = 1:200
    r = sin(c/200*2*pi) * 99
    mat(round(r) + 100,c) = 1;
end

myplot(mat);

[H,Theta,Rho] = hough(mat,'RhoResolution',1);

myplot(Theta,Rho,H); 
hold on;
peaks = houghpeaks(H, 3, 'Threshold', 0.3 * max(H(:)));
for k = 1:size(peaks,1)
   xy = [Theta(peaks(k,2)), Rho(peaks(k,1))];
   plot3(xy(1),xy(2),max(H(:)),'s','LineWidth',2);
end