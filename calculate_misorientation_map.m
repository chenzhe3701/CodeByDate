% calculate misorientation, using phi1, phi, phi2 maps.
% inds indicate the points to do the calculation
% euler_ref is the reference orientation 1x3
% 2021-10-17

function map = calculate_misorientation_map(phi1, phi, phi2, euler_ref, inds)
    [nR,nC] = size(phi1);
    map = zeros(size(phi1)) * nan;
    for iR = 1:nR
       for iC = 1:nC 
           if inds(iR,iC) == 1
               euler = [phi1(iR,iC), phi(iR,iC), phi2(iR,iC)];
               map(iR,iC) = calculate_misorientation_euler_d(euler, euler_ref, 'hcp');
           end
       end
    end
end