function [phi] = kernel(r)
% KERNEL radial basis kernel
%
% INPUTS:
% <r> = euclidean distance between parameter locations
%
% OUTPUTS:
% <phi> = kernel value

% thin plate spline
phi = (r^2) * log(r);

end
