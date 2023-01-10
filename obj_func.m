function [J] = obj_func(w,dx)
% OBJ_FUNC discretized version of the functional (or objective function) J
%
% INPUTS:
% <w> = steady-state solution
% <dx> = length of discretization grid
%
% OUTPUTS:
% <J> = functional value for steady-state solution w

J = sum(w*dx);

end

