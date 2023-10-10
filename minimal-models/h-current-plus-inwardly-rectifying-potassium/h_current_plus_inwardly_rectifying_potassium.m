function integrand = h_current_plus_inwardly_rectifying_potassium(X, I, C, gL, EL, gKir, EK, gh, Eh, VhKir, khKir, Vh, kh, C_base, C_amp, V_max, sigma)
% h_current_plus_inwardly_rectifying_potassium create a function handle of h-current plus inwardly rectifying potassium model.
%
% integrand = h_current_plus_inwardly_rectifying_potassium(X, I, C, gL, EL, gKir, EK, gh, Eh, VhKir, khKir, Vh, kh, C_base, C_amp, V_max, sigma)
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
% C_base, C_amp, V_max, sigma : numeric
%   parameters of voltage-sensitive time constant
%   tau_h = C_base + C_amp.*exp(-((V_max-V)./sigma).^2)
%
% Returns
% -------
% integrand : vector(numeric)
%   h-current plus inwardly rectifying potassium model
%
    hKir_inf = 1 ./ (1 + exp((VhKir-X(1))./khKir));
    h_inf = 1 ./ (1 + exp((Vh-X(1))./kh));
    tau_h = C_base + C_amp.*exp(-((V_max-X(1))./sigma).^2);

    integrand = zeros(2,1);
    integrand(1) = (I - gL*(X(1)-EL) - gKir*hKir_inf*(X(1)-EK) - gh*X(2)*(X(1)-Eh)) / C;
    integrand(2) = (h_inf-X(2)) / tau_h;
end