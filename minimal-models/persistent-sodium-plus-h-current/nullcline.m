function [V_null, h_null] = nullcline(V, I, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh)
% nullcline calculate nullclines.
% 
% [V_null, h_null] = nullcline(V, I, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh)
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
% gh : numeric
%   conductance of h-current [nS]
% Eh : numeric
%   equilibrium potential of h-current [mV]
% Vm, Vh : numeric
% km, kh : numeric
%   parameters of steady-state activation (or inactivation) curves
%   p_inf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
%
% Returns
% -------
% V_null : vector(numeric)
%   V-nullcline
% h_null : vector(numeric)
%   h-nullcline
%
    m_inf = 1 ./ (1 + exp((Vm-V)./km));
    h_inf = 1 ./ (1 + exp((Vh-V)./kh));

    V_null = (I - gL.*(V-EL) - gNa.*m_inf.*(V-ENa)) ./ (gh.*(V-Eh));
    h_null = h_inf;
end