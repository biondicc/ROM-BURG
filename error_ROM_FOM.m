function [J_error] = error_ROM_FOM(w_bar, x, b, ic)
% ERROR_ROM_FOM estimates the error in the functional for solutions from
% the FOM and the ROM
%
% INPUTS:
% <w_bar> = approximate solution from the LSPG projection
% <x> = vector of x-values
% <b> = parameter values
% <ic> = initial condition for the PDE
%
% OUTPUTS:
% <J_error> = error in the functional value between the FOM and the ROM

dx = x(2) - x(1);
[R, dR_dw] = upwind_res(w_bar, x, b, ic);
LHS = dR_dw.';
RHS = ones(length(x), 1)*dx;
adj = LHS \ RHS;
J_error = - (adj.') * R;
    
end

