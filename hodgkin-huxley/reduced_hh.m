function integrand = reduced_hh(X, C, I, gL, EL, gNa, ENa, gK, EK, p)
% reduced_hh create a function handle of reduced Hodgkin-Huxley model.
% 
% integrand = reduced_hh(X, C, I, gL, EL, gNa, ENa, gK, EK, p)
% 
% Parameters
% ----------
% X : vector(numeric)
%   X = [V, m, h, n]
%   V : membrane potential
%   m, h, n : gating variable 
% C : numeric
%   membrane capacitance [μF/cm^2]
% I : numeric
%   external stimulus [pA]
% gL : numeric
%   leakage conductance [nS]
% EL : numeric
%   resting potential [mV]
% gNa : numeric
%   sodium conductance [nS]
% ENa : numeric
%   sodium equilibrium potential [mV]
% gK : numeric
%   potassium conductance [nS]
% EK : numeric
%   potassium equilibrium potential [mV]
% p : vector(numeric)
%   parameters of approximation straight line
%
% Returns
% -------
% integrand : vector(numeric)
%   reduced Hodgkin-Huxley equation
%
    integrand = zeros(4,1);
    [alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n] = gating_variable(X(1));

    tau_m = 1/(alpha_m + beta_m);  m_inf = alpha_m*tau_m;
    tau_h = 1/(alpha_h + beta_h);  h_inf = alpha_h*tau_h;
    tau_n = 1/(alpha_n + beta_n);  n_inf = alpha_n*tau_n;

    integrand(1) = (I - gL*(X(1) - EL) - gNa*(m_inf^3)*(p(2) + p(1)*X(4))*(X(1) - ENa) - gK*(X(4)^4)*(X(1) - EK)) / C;
    integrand(2) = (m_inf - X(2))/tau_m;
    integrand(3) = (h_inf - X(3))/tau_h;
    integrand(4) = (n_inf - X(4))/tau_n;
end