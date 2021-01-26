function [q_out, val_min] = find_closest_quaternion_hcp(q_input, q_ref)
% input 2 quaternions, assume HCP symmetry
% Use the q_ref as reference, find among the crystallographically
% equivalent quaternions of q_input, that is most similar to q_ref without
% considering symmetry.
% 
% chenzhe, 2021-01-22

% HCP symmetry operation.  Derived from code by varying EulerAngles.  To 
% understand, combine with the cubic_symmetry code, and convert to
% angle-axis pairs for better understanding, and better visualization.
S = [1, 0, 0, 0;
    sqrt(3)/2, 0, 0, 1/2;
    1/2, 0, 0, sqrt(3)/2;
    0, 0, 0, -1;
    1/2, 0, 0, -sqrt(3)/2;
    sqrt(3)/2, 0, 0, -1/2;
    0, 1, 0, 0;
    0, sqrt(3)/2, -1/2, 0;
    0, 1/2, -sqrt(3)/2, 0;
    0, 0, 1, 0;
    0, 1/2, sqrt(3)/2, 0;
    0, sqrt(3)/2, 1/2, 0];

% crystallographically equivalent ones
qs = quatmultiply(q_input, S);

% calculate difference
clear thetad;
for ii = 1:size(qs,1)
    thetad(ii,:) = quat2axang(quatmultiply(qs(ii,:), quatconj(q_ref))); 
end
thetad = thetad(:,4)/pi*180;
[val_min, ind_S] = min(abs(thetad));
q_out = qs(ind_S,:);

if  q_out(1) * q_ref(1) < 0
    q_out = -q_out;
end

end



