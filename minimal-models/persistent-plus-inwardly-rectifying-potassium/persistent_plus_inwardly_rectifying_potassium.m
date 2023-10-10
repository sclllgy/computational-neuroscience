function integrand = persistent_plus_inwardly_rectifying_potassium(X, I, C, gKir, EK, gK, Vh, kh, Vn, kn, tau_n)
% persistent_plus_inwardly_rectifying_potassium create a function handle of persistent plus inwardly rectifying potassium model.
% 
% integrand = persistent_plus_inwardly_rectifying_potassium(X, I, C, gKir, EK, gK, Vh, kh, Vn, kn, tau_n)
% 
% Parameters
% ----------
% X : vector(numeric)
%   X = [V, n]
%   V : membrane potential [mV]
%   n : K^+ activation variable
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
% integrand : vector(numeric)
%   persistent plus inwardly rectifying potassium model
%
    h_inf = 1 ./ (1 + exp((Vh-X(1))./kh));
    n_inf = 1 ./ (1 + exp((Vn-X(1))./kn));

    integrand = zeros(2,1);
    integrand(1) = (I - gKir*h_inf*(X(1)-EK) - gK*X(2)*(X(1)-EK)) / C;
    integrand(2) = (n_inf-X(2)) / tau_n;
end