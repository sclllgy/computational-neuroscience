function [alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n] = gating_variable(V)
% gating_variable calculate gating variables, m, h and n.
% 
% [alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n] = gating_variable(V)
% 
% Parameters
% ----------
% V : vector(numeric)
%   membrane potential [mV]
%
% Returns
% -------
% alpha_m, beta_m : vector(numeric)
%   rate constant of sodium channel
% alpha_h, beta_h : vector(numeric)
%   rate constant of sodium channel
% alpha_n, beta_n : vector(numeric)
%   rate constant of potassium channel
%
    alpha_m = 0.1*(25 - V) ./ (exp((25-V)/10) - 1);
    beta_m  = 4*exp(-V/18);
    alpha_h = 0.07*exp(-V/20);
    beta_h  = 1 ./ (exp((30-V)/10) + 1);
    alpha_n = 0.01.*(10 - V) ./ (exp((10-V)/10) - 1);
    beta_n  = 0.125.*exp(-V/80);
end