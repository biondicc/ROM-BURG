function [exp_s] = source(x, mu)
% SOURCE source term in the PDE and residual equation
%
% INPUTS:
% <x> = vector of x-values
% <mu> = parameter values
%
% OUTPUTS:
% <exp_s> = vector of source function evaluated a x-values

if length(mu) > 1
    b = mu(1);
    a = mu(2);
else
    b = mu;
    a = 0.02;
end

exp_s = a*exp(b*x);
%exp_s = exp(b*x) + a*(x/10)^3 - 0.01*x^2;

end

