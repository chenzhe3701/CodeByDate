function [SIR,tVec] = ssir_model(tspan,X0)

% initla S,I,R,t
t = 0;    
tVec = t;

S = X0(1);
I = X0(2);
R = X0(3);
SIR = [S,I,R];

nCount = 0;
while t < tspan
    nCount = nCount + 1;
    gamma_t = 0.0247 * exp(-0.111*t);
    
    beta0 = 7.58e-7;
    beta_t = beta0 * exp(-0.0506*t);
    
    % step 1, generate random numbers
    r_1 = rand();
    r_2 = rand();
    
    % step 2, compute alpha values
    alpha_1 = S * I * beta_t;
    alpha_2 = I * gamma_t;
    
    % step 3, alpha_0
    alpha_0 = alpha_1 + alpha_2;
    
    % step 4, tau
    tau = 1/alpha_0 * log(1/r_1);
    
    % step 5, compute S,I,R at t+tau, and update S,I,R
    if r_2 < alpha_1/alpha_0
        S = S - 1;
        I = I + 1;
        R = R;
    elseif r_2 >= alpha_1/alpha_0
        S = S;
        I = I - 1;
        R = R + 1;
    end
    
    % step 6, check
    if I == 0
        I = 1;
        S = S-1;
    end
    
    t = tVec(end) + tau;
    tVec = [tVec; t];
    SIR = [SIR; S,I,R];
    
end

end