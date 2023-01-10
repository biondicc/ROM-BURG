function [error] = error_from_RBF(mu_new, mu_orig, w)
% ERROR_FROM_RBF estimate error between FOM and ROM at new parameter 
% location using the Radial Basis Function
%
% INPUTS:
% <mu_new> = vector containing design parameters
% <mu_original> = locations of design parameters used to build the RBF
% <w> = weights for RBF
%
% OUTPUTS:
% <error> = estimated error between FOM and ROM at new parameter location

phi = zeros(1,length(mu_orig));
for i=1:length(mu_orig)
    mu_i = mu_orig(i,:);
    r = norm((mu_i-mu_new),2);
    if r == 0
        phi(1,i) = 0;
    else
        phi(1,i) = kernel(r);
    end
end

error = phi * w;

end

