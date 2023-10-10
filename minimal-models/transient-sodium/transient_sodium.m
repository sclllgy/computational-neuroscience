function integrand = transient_sodium(X, I, C, gL, EL, gNa, ENa, Vm, km, Vh, kh, tau_h)
% transient_sodium create a function handle of transient sodium model.
% 
% integrand = transient_sodium(X, I, C, gL, EL, gNa, ENa, Vm, km, Vh, kh, tau_h)
% 
% Parameters
% ----------
% X : vector(numeric)
%   X = [V, h]
%   V : membrane potential [mV]
%   h : Na^+ inactivation variable
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
% integrand : vector(numeric)
%   transient sodium model
%
    m_inf = 1 ./ (1 + exp((Vm-X(1))./km));
    h_inf = 1 ./ (1 + exp((Vh-X(1))./kh));

    integrand = zeros(2,1);
    integrand(1) = (I - gL*(X(1)-EL) - gNa*(m_inf^3)*X(2)*(X(1)-ENa)) / C;
    integrand(2) = (h_inf-X(2)) / tau_h;
end