function [V_null, n_null] = nullcline(V, I, gKir, EK, gK, Vh, kh, Vn, kn)
% nullcline calculate nullclines.
% 
% [V_null, h_null] = nullcline(V, I, gKir, EK, gK, Vh, kh, Vn, kn)
% 
% Parameters
% ----------
% V : vector(numeric)
%   membrane potential [mV]
% I : numeric
%   external stimulus [pA]
% C : nureric
%   membrane capacitance [μF]
% gKir : numeric
%   inwardly rectifying potassium conductance [nS]
% EK : numeric
%   sodium equilibrium potential [mV]
% gK : numeric
%   sodium conductance [nS]
% Vh, Vn : numeric
% kh, kn : numeric
%   parameters of steady-state activation (or inactivation) curves
%   p_inf = 1./ (1 + (exp(Vp-V)./kp)), p = h or n
%
% Returns
% -------
% V_null : vector(numeric)
%   V-nullcline
% n_null : vector(numeric)
%   n-nullcline
%
    h_inf = 1 ./ (1 + exp((Vh-V)./kh));
    n_inf = 1 ./ (1 + exp((Vn-V)./kn));

    V_null = I ./ (gK.*(V-EK)) - gKir.*h_inf ./ gK;
    n_null = n_inf;
end