function [V_null, n_null] = nullcline(V, I, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn)
% nullcline calculate nullclines.
% 
% [V_null, n_null] = nullcline(V, I, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn)
% 
% Parameters
% ----------
% V : vector(numeric)
%   membrane potential [mV]
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
%
% Returns
% -------
% V_null : vector(numeric)
%   V-nullcline
% n_null : vector(numeric)
%   n-nullcline
%
    m_inf = 1 ./ (1 + exp((Vm-V)./km));
    n_inf = 1 ./ (1 + exp((Vn-V)./kn));

    V_null = (I - gL.*(V-EL) - gNa.*m_inf.*(V-ENa)) ./ (gK.*(V-EK));
    n_null = n_inf;
end