function integrand = persistent_sodium_plus_potassium(X, I, C, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn, tau_n)
% persistent_sodium_plus_potassium create a function handle of persistent sodium plus potassium model.
% 
% integrand = persistent_sodium_plus_potassium(X, I, C, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn, tau_n)
% 
% Parameters
% ----------
% X : vector(numeric)
%   X = [V, n]
%   V : membrane potential [mV]
%   n : K^+ activation variable
% I : numeric
%   external stimulus [pA]
% C : nureric
%   membrane capacitance [μF]
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
% Vm, Vn : numeric
% km, kn : numeric
%   parameters of steady-state activation (or inactivation) curves
%   p_inf = 1./ (1 + (exp(Vp-V)./kp)), p = m or n
% tau_n : numeric
%   time constant of n_inf [ms]
%
% Returns
% -------
% integrand : vector(numeric)
%   persistent sodium plus potassium model
%
    m_inf = 1 ./ (1 + exp((Vm-X(1))./km));
    n_inf = 1 ./ (1 + exp((Vn-X(1))./kn));

    integrand = zeros(2,1);
    integrand(1) = (I - gL*(X(1)-EL) - gNa*m_inf*(X(1)-ENa) - gK*X(2)*(X(1)-EK)) / C;
    integrand(2) = (n_inf - X(2)) / tau_n;
end