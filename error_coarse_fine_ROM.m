function [J_error] = error_coarse_fine_ROM(w_bar, w_bar_h, x, b, ic, V_h)
%ERROR_COARSE_FINE_ROM Summary of this function goes here
%   Detailed explanation goes here
dx = x(2) - x(1);
[R, dR_dw] = upwind_res(w_bar, x, b, ic);
[~, dR_dw_h] = upwind_res(w_bar_h, x, b, ic);
LHS = (V_h.') * (dR_dw.') * dR_dw * V_h;
LHS = LHS.';
dJ_dw = ones(length(x), 1)*dx;
RHS = dJ_dw.' * V_h;
RHS = RHS.';
adj = LHS \ RHS;
W_T = (V_h.') * (dR_dw_h.');
R_h = (W_T) *  R;
J_error = - (adj.') * R_h;
end

