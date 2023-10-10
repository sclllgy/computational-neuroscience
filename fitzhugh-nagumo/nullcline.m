function [V_null, w_null] = nullcline(V, I, a, b, c)
% nullcline calculate nullclines.
% 
% [V_null, w_null] = nullcline(V, I, a, b, c)
% 
% Parameters
% ----------
% V : vector(numeric)
%   membrane potential [mV]
% I : numeric
%   external stimulus [pA]
% a : numeric
% b : numeric
% c : numeric
% 
% Returns
% -------
% V_null : vector(numeric)
%   V-nullcline
% w_null : vector(numeric)
%   w-nullcline
%
    V_null = V.*(a - V).*(V - 1) + I;
    w_null = (b/c).*V;
end