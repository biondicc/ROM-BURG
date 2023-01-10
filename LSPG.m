function [w_bar, iter] = LSPG(V, w_ref, tol, x, b, ic)
% LSPG use least-squares Petrov-Galerkin projection to find approximate
% solution using the POD basis
%
% INPUTS:
% <V> = N by M basis matrix from POD
% <w_ref> = reference solution (average of snapshots)
% <tol> = error tolerance
% <x> = vector of x-values
% <b> = parameter values
% <ic> = initial condition for the PDE
%
% OUTPUTS:
% <w_bar> = approximate solution at specified parameter location from ROM
% <iter> = number of Newton iterations to converge

tol = 10^-6;

w_0 = ones(length(V),1);
w_bar = w_ref + V * (V.' * (w_0 - w_ref));

[R, dR_dw] = upwind_res(w_bar, x, b, ic);
r = norm(V.' * dR_dw.' * R,2);
a = 0.05;
iter = 1;
while r > tol
    [R, dR_dw] = upwind_res(w_bar, x, b, ic);
    RHS = - (V.') * (dR_dw.') * R;
    LHS = (V.') * (dR_dw.') * dR_dw * V;
    p_k = LHS \ RHS;
    w_bar = w_bar + V * (a * p_k);
    r = norm(V.' * dR_dw.' * R,2);
    iter = iter + 1;
end

end

