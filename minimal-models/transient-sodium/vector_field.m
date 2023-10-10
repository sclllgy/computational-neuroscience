function [dVdt, dhdt] = vector_field(V, h, I, C, gL, EL, gNa, ENa, Vm, km, Vh, kh, tau_h)
% vector_field calculate vector field.
% 
% [dVdt, dhdt] = calc_vector_field(V, h, I, C, gL, EL, gNa, ENa, Vm, km, Vh, kh, tau_h)
% 
% Parameters
% ----------
% V : vector(numeric)
%   membrane potential [mV]
% h : vector(numeric)
%   Na^+ inactivation variable
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
% Vm, Vh : numeric
% km, kh : numeric
%   parameters of steady-state activation (or inactivation) curves
%   p_inf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
% tau_h : numeric
%   time constant of h_inf [ms]
%
% Returns
% -------
% dVdt : vector(numeric)
%   time derivative of V
% dhdt : vector(numeric)
%   time derivative of h
%
    m_inf = 1 ./ (1 + exp((Vm-V)./km));
    h_inf = 1 ./ (1 + exp((Vh-V)./kh));

    dVdt = (I - gL.*(V-EL) - gNa.*(m_inf.^3).*h.*(V-ENa)) ./ C;
    dhdt = (h_inf-h) ./ tau_h;
end