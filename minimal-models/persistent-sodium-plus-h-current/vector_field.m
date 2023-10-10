function [dVdt, dhdt] = vector_field(V, h, I, C, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh, C_base, C_amp, V_max, sigma)
% vector_field calculate vector field.
% 
% [dVdt, dhdt] = vector_field(V, h, I, C, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh, C_base, C_amp, V_max, sigma)
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
% gh : numeric
%   conductance of Ih [nS]
% Eh : numeric
%   equilibrium potential of Ih [mV]
% Vm, Vh : numeric
% km, kh : numeric
%   parameters of steady-state activation (or inactivation) curves
%   p_inf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
% C_base, C_amp, V_max, sigma : numeric
%   parameters of voltage-sensitive time constant [ms]
%   tau_h = C_base + C_amp.*exp(-((V_max-V)./sigma).^2)
%
% Returns
% -------
% dVdt : vector(numeric)
%   derivative function of V
% dhdt : vector(numeric)
%   derivative function of h
%
    m_inf = 1 ./ (1 + exp((Vm-V)./km));
    h_inf = 1 ./ (1 + exp((Vh-V)./kh));
    tau_h = C_base + C_amp.*exp(-((V_max-V)./sigma).^2);

    dVdt = (I - gL.*(V-EL) - gNa.*m_inf.*(V-ENa) - gh.*h.*(V-Eh)) ./ C;
    dhdt = (h_inf-h) ./ tau_h;
end