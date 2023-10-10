function integrand = fitzhugh_nagumo(X, I, a, b, c)
% fitzhugh_nagumo create a function handle of FitzHugh-Nagumo model.
% 
% integrand = fitzhugh_nagumo(X, I, a, b, c)
% 
% Parameters
% ----------
% X : vector(numeric)
%   X = [dVdt, dwdt]
%   dVdt : time derivative of V
%   dwdt : time derivative of w
% I : numeric
%   external stimulus [pA]
%
% Returns
% -------
% integrand : vector(numeric)
%   FitzHugh-Nagumo model
%
    integrand = zeros(2,1);
    integrand(1) = X(1).*(a - X(1)).*(X(1) - 1) - X(2) + I;
    integrand(2) = b.*X(1) - c.*X(2);
end
