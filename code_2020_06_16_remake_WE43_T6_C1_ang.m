% remake .ang file for WE43_T6_C1, so the phase data becomes readable for
% Mtex
% 2020-06-16

fName = 'D:\WE43_T6_C1\EBSD Data\WE43_T6_C1_Stitched_Dialated.ang';
[phi1,phi,phi2,x,y,IQ,CI,Phase,Intensity,Fit] = read_ang(fName);

Phase = Phase + 1;
Phase(1:2,1:2)=0;   % special treatment so that Mtex can read
phi1(1) = 0.2; % rad, ~= 11.46 deg
phi(1) = 0.1;   % 5.73 deg
phi2(1) = 0.05; % 2.87 deg

data.phi1 = phi1;
data.phi = phi;
data.phi2 = phi2;
data.x = x;
data.y = y;
data.IQ = IQ;
data.CI = CI;
data.Phase = Phase;
data.Intensity = Intensity;
data.Fit = Fit;

fName_target = 'D:\WE43_T6_C1\EBSD Data\WE43_T6_C1_Stitched_Dialated_rewrite.ang';
write_ang(fName_target, 'Mg', data)