function [V, w_ref] = compute_POD_basis(S)
% COMPUTE_POD_BASIS find the POD basis from solution snapshots
%
% INPUTS:
% <S> = N by M matrix of solution snapshots
%
% OUTPUTS:
% <V> = N by M basis matrix from "thin" SVD
% <w_ref> = reference solution (average of snapshots)

w_ref = mean(S, 2);
Q = S - w_ref;
[V,~,~] = svd(Q,"econ");

end

