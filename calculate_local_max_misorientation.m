% calculate the max misorientation for each pixel with its four adjacent pixels, for HCP
% 2021-10-17
function misorientation_max = calculate_local_max_misorientation(phi1_local, phi_local, phi2_local, CS)

% crystal symmetry
if ~exist('CS','var')
    CS = 'HCP';
end

[nR,nC] = size(phi1_local);
misorientation_max = zeros(nR,nC);
for iR = 1:nR
    for iC = 1:nC
        neighbor_orientation = [];
        % add upper neighbor, bottom neighbor, left neighbor, right neighbor
        if iR>1
            neighbor_orientation = [neighbor_orientation; phi1_local(iR-1,iC), phi_local(iR-1,iC), phi2_local(iR-1,iC)];
        end
        if iR<nR
            neighbor_orientation = [neighbor_orientation; phi1_local(iR+1,iC), phi_local(iR+1,iC), phi2_local(iR+1,iC)];
        end
        if iC>1
            neighbor_orientation = [neighbor_orientation; phi1_local(iR,iC-1), phi_local(iR,iC-1), phi2_local(iR,iC-1)];
        end
        if iC<nC
            neighbor_orientation = [neighbor_orientation; phi1_local(iR,iC+1), phi_local(iR,iC+1), phi2_local(iR,iC+1)];
        end
        
        euler_current = [phi1_local(iR,iC), phi_local(iR,iC), phi2_local(iR,iC)];
        
        misorientation = [];
        for ii = 1:size(neighbor_orientation,1)
            misorientation(ii) = calculate_misorientation_euler_d(euler_current, neighbor_orientation(ii,:), CS);
        end
        misorientation_max(iR,iC) = max(misorientation);
    end
end

end