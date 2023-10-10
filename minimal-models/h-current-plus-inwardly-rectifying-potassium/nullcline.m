function [V_null, h_null] = nullcline(V, I, gL, EL, gKir, EK, gh, Eh, VhKir, khKir, Vh, kh)
% nullcline calculate nullclines.
% 
% [V_null, h_null] = nullcline(V, I, gL, EL, gKir, EK, gh, Eh, VhKir, khKir, Vh, kh)
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
% gKir : numeric
%   inwardly rectifying potassium conductance [nS]
% EK : numeric
%   sodium equilibrium potential [mV]
% gh : numeric
%   conductance of h-current [nS]
% Eh : numeric
%   equilibrium potential of h-current [mV]
% VhKir, Vh : numeric
% khKir, kh : numeric
%   parameters of steady-state activation (or inactivation) curves
%   p_inf = 1./ (1 + (exp(Vp-V)./kp)), p = hKir or h
% 
% Returns
% -------
% V_null : vector(numeric)
%   V-nullcline
% h_null : vector(numeric)
%   h-nullcline
%
    hKir_inf = 1 ./ (1 + exp((VhKir-V)./khKir));
    h_inf    = 1 ./ (1 + exp((Vh-V)./kh));

    V_null = (I - gL.*(V-EL) - gKir.*hKir_inf.*(V-EK)) ./ (gh.*(V-Eh));
    h_null = h_inf;
end