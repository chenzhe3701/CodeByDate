function dydt = sir_model(t,X)
% returns the rhs of ODE system of the SIR model

global beta;
gamma_t = 0.0247 * exp(-0.111*t);
beta_t = beta * exp(-0.0506*t);

S = X(1);
I = X(2);
R = X(3);

dSdt = -beta_t * S * I;
dIdt = beta_t * S * I - gamma_t * I;
dRdt = gamma_t * I;

dydt = [dSdt; dIdt; dRdt];

end