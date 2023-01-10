function [mu_max, max_error] = max_error_param(S_mu, P_mu, errors)
% MAX_ERROR_PARAM estimate the parameter location where the error in the
% functional is maximum using radial basis function interpolation
%
% INPUTS:
% <S_mu> = matrix of snapshot parameter locations used to build ROM
% <P_mu> = parameter locations for error estimates/building RBF
% <errors> = error in the functional at each parameter location in P_mu
%
% OUTPUTS:
% <mu_max> = estimated parameter location of maximum error
% <max_error> = maximum error over the parameter field 

A = vertcat(S_mu, P_mu);
B = vertcat(zeros(length(S_mu),1),errors);
A_scaled = zeros(size(A));
A_size = size(A);
num_snaps = A_size(1);
num_params = A_size(2);
for i=1:num_params
    mu_min = min(A(:,i));
    mu_max = max(A(:,i));
    for j = 1:num_snaps
        A_scaled(j,i) = (A(j, i) - mu_min)/ (mu_max - mu_min);
    end
end

[w] = RBF(A_scaled, B);

for k=1:num_params
    mu_min = min(A_scaled(:,k));
    mu_max = max(A_scaled(:,k));
    mu_samples(:,k) = mu_min + (mu_max - mu_min)*rand(1000,1);
end

sample_error = zeros(length(mu_samples),1);
for l=1:length(mu_samples)
    sample_error(l,1) = error_from_RBF(mu_samples(l,:),A_scaled,w);
end
[max_error, mu_idx] = max(abs(sample_error));
mu_max_scaled = mu_samples(mu_idx,:);
for k=1:num_params
    mu_min = min(A(:,k));
    mu_max = max(A(:,k));
    mu_max(1,k) = mu_min + (mu_max - mu_min)*mu_max_scaled(1,k);
end

end

