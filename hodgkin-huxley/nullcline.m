function [V_null, n_null] = nullcline(V, I, gL, EL, gNa, ENa, gK, EK, p)
% nullcline calculate nullclines.
% 
% [V_null, n_null] = nullcline(V, I, gL, EL, gNa, ENa, gK, EK, p)
% 
% Parameters
% ----------
% V : numeric
%   membrane potential [mV]
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
% V_null : numeric
%   V-nullcline
% n_null : numeric
%   n-nullcline
%
    syms n

    [alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n] = gating_variable(V);
    tau_m = 1./(alpha_m + beta_m);  m_inf = alpha_m .* tau_m;
    tau_h = 1./(alpha_h + beta_h);  h_inf = alpha_h .* tau_h;
    tau_n = 1./(alpha_n + beta_n);  n_inf = alpha_n .* tau_n;

    sol_hh = vpasolve(I - gL*(V-EL) - gNa*(m_inf^3)*(p(2)+p(1)*n)*(V-ENa) - gK*(n^4)*(V-EK) == 0, n, [-Inf Inf]);
    V_null = double(max(sol_hh));
    n_null = n_inf;
end