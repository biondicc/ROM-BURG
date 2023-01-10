function [r, drdw] = upwind_res(w, x, b, ic)
% UPWIND_RES 
%
% INPUTS:
% <w> = approximate steady-state solution
% <x> = vector of x-values
% <b> = parameter values
% <ic> = initial condition for the PDE
%
% OUTPUTS:
% <r> = residual at approximate solution
% <drdw> = derivative of the residual

dx = x(2) - x(1);
r = zeros(length(x), 1);
for i=1:length(x)
    if i == 1
        r(i) = (w(i)^2 - ic^2)/(2*dx) - source(x(i),b);
    else
        r(i) = (w(i)^2 - w(i-1)^2)/(2*dx) - source(x(i),b);
    end
end

drdw = zeros(length(x), length(x));
for i=1:length(x)
    if i == 1
        drdw(1,1) = w(i)/dx;
    else
        drdw(i,i-1) = -w(i-1)/dx;
        drdw(i,i) = w(i)/dx;
    end
end
end

