function integrand = persistent_sodium_plus_h_current(X, I, C, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh, C_base, C_amp, V_max, sigma)
% persistent_sodium_plus_h_current create a function handle of persistent sodium plus h-current model.
% 
% integrand = persistent_sodium_plus_h_current(X, I, C, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh, C_base, C_amp, V_max, sigma)
%
% Parameters
% ----------
% X : vector(numeric)
%   X = [dVdt, dhdt]
%   dVdt : time derivative of V
%   dhdt : time derivative of h
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
% C_base, C_amp, V_max, sigma : numeric
%   parameters of voltage-sensitive time constant [ms]
%   tau_h = C_base + C_amp.*exp(-((V_max-V)./sigma).^2)
%
% Returns
% -------
% integrand : vector(numeric)
%   persistent sodium plus h-current model
%
    m_inf = 1 ./ (1 + exp((Vm-X(1))./km));
    h_inf = 1 ./ (1 + exp((Vh-X(1))./kh));
    tau_h = C_base + C_amp.*exp(-((V_max-X(1))./sigma).^2);

    integrand = zeros(2,1);
    integrand(1) = (I - gL*(X(1)-EL) - gNa*m_inf*(X(1)-ENa) - gh*X(2)*(X(1)-Eh)) / C;
    integrand(2) = (h_inf-X(2)) / tau_h;
end