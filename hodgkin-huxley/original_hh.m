function integrand = original_hh(X, C, I, gL, EL, gNa, ENa, gK, EK)
% original_hh create a function handle of original Hodgkin-Huxley model.
% 
% integrand = original_hh(X, C, I, gL, EL, gNa, ENa, gK, EK)
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
%
% Returns
% -------
% integrand : vector(numeric)
%   original Hodgkin-Huxley equation
%
    integrand = zeros(4,1);
    [alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n] = gating_variable(X(1));

    integrand(1) = (I - gL*(X(1) - EL) - gNa*(X(2)^3)*X(3)*(X(1) - ENa) - gK*(X(4)^4)*(X(1) - EK)) / C;
    integrand(2) = alpha_m*(1 - X(2)) - beta_m*X(2);
    integrand(3) = alpha_h*(1 - X(3)) - beta_h*X(3);
    integrand(4) = alpha_n*(1 - X(4)) - beta_n*X(4);
end