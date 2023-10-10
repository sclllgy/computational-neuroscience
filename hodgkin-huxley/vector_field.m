function [dVdt, dndt] = vector_field(V, n, C, I, gL, EL, gNa, ENa, gK, EK, p)
% vector_field calculate vector field.
% 
% [dVdt, dndt] = vector_field(V, n, C, I, gL, EL, gNa, ENa, gK, EK, p)
% 
% Parameters
% ----------
% V : vector(numeric)
%   membrane potential [mV]
% n : vector(numeric)
%   K^+ activation variable
% C : numeric
%   membrane capacitance [μF]
% I : numeric
%   external stimulus [pA]
% gL : numeric
%   leakage conductance [nS]
% EL : numeric
%   resting potential [mV]
% gNa : numeric
%   sodium conductance [nS]
% ENa : numeric
%   sodium equilibrium potential [mV]
% gK : numeric
%   potassium conductance [nS]
% EK : numeric
%   potassium equilibrium potential [mV]
% p : vector(numeric)
%   parameters of approximation straight line
%
% Returns
% -------
% dVdt : vector(numeric)
%   time derivative of V
% dndt : vector(numeric)
%   time derivative of n
%
    [alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n] = gating_variable(V);
    tau_m = 1./(alpha_m + beta_m);  m_inf = alpha_m.*tau_m;
    tau_h = 1./(alpha_h + beta_h);  h_inf = alpha_h.*tau_h;
    tau_n = 1./(alpha_n + beta_n);  n_inf = alpha_n.*tau_n;

    dVdt = (I - gL.*(V - EL) - gNa.*(m_inf.^3).*(p(2) + p(1).*n).*(V - ENa) - gK*(n.^4).*(V - EK)) ./ C;
    dndt = (n_inf - n) ./ tau_n;
end