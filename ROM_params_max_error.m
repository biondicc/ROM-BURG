function [P_mu_star] = ROM_params_max_error(S_params, mu_max_error)
% ROM_PARAMS_MAX_ERROR find the new ROM probe points associated with the
% parameter location of maximum error which is going to be added to the
% snapshots
%
% INPUTS:
% <S_params> = matrix of snapshot parameter locations
% <mu_max_error> = current estimated parameter location of maximum error
%
% OUTPUTS:
% <P_mu_star> = (ROM probes) parameter locations for error estimates/building RBF

S_scaled = zeros(size(S_params));
S_size = size(S_params);
num_snaps = S_size(1);
num_params = S_size(2);
for i=1:num_params
    mu_min = min(S_params(:,i));
    mu_max = max(S_params(:,i));
    for j = 1:num_snaps
        S_scaled(j,i) = (S_params(j, i) - mu_min)/ (mu_max - mu_min);
    end
    mu_max_error_scaled(1,i) = (mu_max_error(1, i) - mu_min)/ (mu_max - mu_min);
end
P_mu_star = [];

for n = 1:num_snaps
    distances(n) = norm((S_scaled(n,:)-mu_max_error_scaled),2);
end
[~,idx] = sort(distances);
n_p = 1;
mu_mid_scaled = [];
mu_mid = [];
while (n_p <= (num_params+1)) && (n_p <= max(idx))
    mu_close = S_scaled(idx(n_p),:);
    mu_mid_scaled(n_p, :) = (mu_close + mu_max_error_scaled)/2;
    n_p = n_p + 1;
end
for i=1:num_params
    mu_min = min(S_params(:,i));
    mu_max = max(S_params(:,i));
    for j = 1:length(mu_mid_scaled)
        mu_mid(j, i) = mu_mid_scaled(j, i)*(mu_max-mu_min) + mu_min;
    end
end
for s = 1:length(mu_mid)
    if isempty(P_mu_star)
        P_mu_star = [P_mu_star, mu_mid(s)];
    elseif num_params == 1
        if ~(ismember(mu_mid(s), P_mu_star)) && (~(ismember(mu_mid(s), S_params)))
           P_mu_star = [P_mu_star, mu_mid(s)];
        end
    else
        if (~(ismember(mu_mid(s), P_mu_star, "rows"))) && (~(ismember(mu_mid(s), S_params, "rows")))
            P_mu_star = [P_mu_star, mu_mid(s)];
        end
    end
end

P_mu_star = P_mu_star.';
        
end

