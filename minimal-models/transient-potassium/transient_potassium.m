function integrand = transient_potassium(X, I, C, gL, EL, gA, EK, Vm, km, Vh, kh, tau_m)
% transient_potassium create a function handle of transient potassium (A-current) model.
% 
% integrand = transient_potassium(X, I, C, gL, EL, gA, EK, Vm, km, Vh, kh, tau_m)
% 
% Parameters
% ----------
% X : vector(numeric)
%   X = [dVdt, dhdt]
%   dVdt : time derivative of V
%   dmdt : time derivative of m
% I : numeric
%   external stimulus [pA]
% C : nureric
%   membrane capacitance [μF]
% gL : numeric
%   leakage conductance [nS]
% EL : numeric
%   resting potential [mV]
% gA : numeric
%   potassium conductance [nS]
% EK : numeric
%   potassium equilibrium potential [mV]
% Vm, Vh : numeric
% km, kh : numeric
%   parameters of steady-state activation (or inactivation) curves
%   p_inf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
% tau_m : numeric
%   time constant of m_inf [ms]
%
% Returns
% -------
% integrand : vector(numeric)
%   transient potassium (A-current) model
%
    m_inf = 1 ./ (1 + exp((Vm-X(1))./km));
    h_inf = 1 ./ (1 + exp((Vh-X(1))./kh));

    integrand = zeros(2,1);
    integrand(1) = (I - gL*(X(1)-EL) - gA*X(2)*h_inf*(X(1)-EK)) / C;
    integrand(2) = (m_inf-X(2)) / tau_m;
end