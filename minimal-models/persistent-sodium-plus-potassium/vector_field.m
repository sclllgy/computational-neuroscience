function [dVdt, dndt] = vector_field(V, n, I, C, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn, tau_n)
% vector_field calculate vector field.
% 
% [dVdt, dndt] = vector_field(V, n, I, C, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn, tau_n)
% 
% Parameters
% ----------
% V : vector(numeric)
%   membrane potential [mV]
% n : vector(numeric)
%   K^+ activation variable
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
% V_null : vector(numeric)
%   V-nullcline
% n_null : vector(numeric)
%   n-nullcline
%
    m_inf = 1 ./ (1 + exp((Vm-V)./km));
    n_inf = 1 ./ (1 + exp((Vn-V)./kn));

    dVdt = (I - gL.*(V-EL) - gNa.*m_inf.*(V-ENa) - gK.*n.*(V-EK)) ./ C;
    dndt = (n_inf-n) ./ tau_n;
end