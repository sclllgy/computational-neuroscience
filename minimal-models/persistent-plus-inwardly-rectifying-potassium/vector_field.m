function [dVdt, dndt] = vector_field(V, n, I, C, gKir, EK, gK, Vh, kh, Vn, kn, tau_n)
% vector_field calculate vector field.
% 
% [dVdt, dndt] = vector_field(V, n, I, C, gKir, EK, gK, Vh, kh, Vn, kn, tau_n)
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
% tau_n : numeric
%   time constant of n_inf [ms]
%
% Returns
% -------
% dVdt : vector(numeric)
%   time derivative of V
% dndt : vector(numeric)
%   time derivative of n
%
    h_inf = 1 ./ (1 + exp((Vh-V)./kh));
    n_inf = 1 ./ (1 + exp((Vn-V)./kn));

    dVdt = (I - gKir.*h_inf.*(V-EK) - gK.*n.*(V-EK)) ./ C;
    dndt = (n_inf-n) ./ tau_n;
end