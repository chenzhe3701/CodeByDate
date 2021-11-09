% In the paper, Table I seems to be problematic.
% This code generates orientations similar to those in the paper, and look
% at what the Schmid factor values should be.

close all;

eulers_d = [0 10 15];   % {10-12} specimen
eulers_d = [0 30 25];   % {11-22} specimen
eulers_d = [0 50 5];    % {0001} specimen
eulers_d = [0 80 2];    % {-1010} specimen
eulers_d = [0 85 28];   % {1-100} specimen

phi_sys = [0 0 0];
phi_error = [0 0 0];
sample_direction_of_interest = [0 0 1];
ss_of_interest = 1;
stressTensor = [0 0 0; 0 0 0; 0 0 -1];
str_material = 'Ti';
str_twin = 'ctwin';

plot_on_IPF(eulers_d, phi_sys, phi_error, sample_direction_of_interest, ss_of_interest, stressTensor, str_material, str_twin)