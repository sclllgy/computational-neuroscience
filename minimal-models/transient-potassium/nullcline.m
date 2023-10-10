function [V_null, m_null] = nullcline(V, I, gL, EL, gA, EK, Vm, km, Vh, kh)
% nullcline calculate nullclines.
% 
% [V_null, m_null] = nullcline(V, I, gL, EL, gA, EK, Vm, km, Vh, kh)
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
% gA : numeric
%   potassium conductance [nS]
% EK : numeric
%   potassium equilibrium potential [mV]
% Vm, Vh : numeric
% km, kh : numeric
%   parameters of steady-state activation (or inactivation) curves
%   p_inf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
%
% Returns
% -------
% V_null : vector(numeric)
%   V-nullcline
% m_null : vector(numeric)
%   m-nullcline
%
    m_inf = 1 ./ (1 + exp((Vm-V)./km));
    h_inf = 1 ./ (1 + exp((Vh-V)./kh));

    V_null = (I - gL.*(V-EL)) ./ (gA.*h_inf.*(V-EK));
    m_null = m_inf;
end