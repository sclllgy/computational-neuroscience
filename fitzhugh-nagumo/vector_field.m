function [dVdt, dwdt] = vector_field(V, w, I, a, b, c)
% vector_field calculate vector field.
% 
% [dVdt, dwdt] = nullcline(V, I, a, b, c)
% 
% Parameters
% ----------
% V : vector(numeric)
%   membrane potential [mV]
% w : vector(numeric)
%   recovery variable
% I : numeric
%   external stimulus [pA]
% a : numeric
% b : numeric
% c : numeric
% 
% Returns
% -------
% dVdt : vector(numeric)
%   time derivative of V
% dwdt : vector(numeric)
%   time derivative of w
%
    dVdt = V.*(a - V).*(V - 1) - w + I;
    dwdt = b.*V - c.*w;
end