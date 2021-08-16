% study the local Schmid Factor metric
% in paper: In situ investigation of extension twinning-detwinning and its effect on the mechanical behavior of AZ31B magnesium alloy  
% B.L.Wu et al, Materials and Design 132(2017)57-65 

close all;
stress = [1 0 0; 0 0 0; 0 0 0];
for ii=1:1000
    % euler of grain 1
    euler_1 = rand(1,3).*[360,180,360];
    euler_2 = rand(1,3).*[360,180,360];
    
    M1 = angle2dcm(euler_1(1)/180*pi, euler_1(2)/180*pi, euler_1(3)/180*pi, 'zxz');
    M2 = angle2dcm(euler_2(1)/180*pi, euler_2(2)/180*pi, euler_2(3)/180*pi, 'zxz');
    
    % Define 'Local Schmid Factor' Fl = cos(lambda) sqrt(1-(cos(alpha))^2)
    % lambda = angle <eta0 = twin shear dir g1, s0 = accommodative shear dir g2>
    % alpha = angle <eta0, n0 = accommodative shear plane normal g2>
    sss = define_SS_cart('Mg','twin');
    
    iss_g1 = 19; %   twin system in g1, should be member of 19:24 to be twin
    iss_g2 = 1; %   twin system in g2,
    
    n1 = sss(1,:,iss_g1);  % g1 plane normal in crystal coordinate
    n1 = n1 * M1;       % transform to lab coordinate
    b1 = sss(2,:,iss_g1);  % g1 direction in crystal coordinate
    b1 = b1 * M1;       % transform to lab coordinate
    
    
    n2 = sss(1,:,iss_g2);  % g1 plane normal in crystal coordinate
    n2 = n2 * M2;       % transform to lab coordinate
    b2 = sss(2,:,iss_g2);  % g1 direction in crystal coordinate
    b2 = b2 * M2;       % transform to lab coordinate
    
    % if g2 is slip, first part always positive
    if iss_g2 <= 18
        F_local = abs(dot(b1, b2)) * sqrt(1-(dot(b1,n2))^2);
    else
        F_local = dot(b1, b2) * sqrt(1-(dot(b1,n2))^2);
    end
        
    F_twin = n1 * stress * b1(:);   % twin Schmid factor
    
    % A < stress/tau_ct  should be satisfied if deformation in g2 can% accommodate twin in g1  
    % tau_cp: CRSS of accommodating system
    % tau_ct: CRSS of extension twin
    tau_ct = 2.4;
    if ismember(iss_g2, 1:3)
        tau_cp = 0.64;
    elseif ismember(iss_g2, 4:6)
        tau_cp = 39.2;
    elseif ismember(iss_g2, 19:24)
        tau_cp = 2.4;
    end
    A = 1/F_twin/F_local * tau_cp/tau_ct + 1/F_twin;  
    
    As(ii) = A;
    F_locals(ii) = F_local;
end

figure; 
histogram(As, [-inf,-100:10:100,inf]); 
title('A distribution');

figure; 
histogram(F_locals); 
title('F-local distribution');




