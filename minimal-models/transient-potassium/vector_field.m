function [dVdt, dmdt] = vector_field(V, m, I, C, gL, EL, gA, EK, Vm, km, Vh, kh, tau_m)
% vector_field calculate vector field.
% 
% [dVdt, dmdt] = vector_field(V, m, I, C, gL, EL, gA, EK, Vm, km, Vh, kh, tau_m)
% 
% Parameters
% ----------
% V : vector(numeric)
%   membrane potential [mV]
% m : vector(numeric)
%   K^+ activation variable
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
% dVdt : vector(numeric)
%   time derivative of V
% dmdt : vector(numeric)
%   time derivative of m
%
    m_inf = 1 ./ (1 + exp((Vm-V)./km));
    h_inf = 1 ./ (1 + exp((Vh-V)./kh));

    dVdt = (I - gL.*(V-EL) - gA.*m.*h_inf.*(V-EK)) ./ C;
    dmdt = (m_inf-m) ./ tau_m;
end