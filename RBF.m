function [w] = RBF(MU, errors)
% RBF find the weights for the RBF interpolation
%
% INPUTS:
% <MU> = matrix of all parameter locations to be used in the RBF
% <errors> = corresponding error estimates at each location in MU
%
% OUTPUTS:
% <w> = weights for the RBF interpolation

A = zeros(length(MU),length(MU));
for i=1:length(MU)
    mu_i = MU(i,:);
    for j =1:length(MU)
        mu_j = MU(j,:);
        r = norm((mu_i-mu_j),2);
        if r == 0
            A(i,j) = 0;
        else
            A(i,j) = kernel(r);
        end
    end
end

w = A \ errors;

end

