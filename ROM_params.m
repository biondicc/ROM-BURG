function [P_mu] = ROM_params(S_params)
% ROM_PARAMS compute ROM probe locations (mid-points) from snapshot
% locations used to build the POD basis
%
% INPUTS:
% <S_params> = matrix of snapshot parameter locations
%
% OUTPUTS:
% <P_mu> = (ROM probes) parameter locations for error estimates/building RBF

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
end
P_mu = [];
for m = 1:num_snaps
    mu_s = S_scaled(m,:);
    for n = 1:num_snaps
        distances(n) = norm((S_scaled(n,:)-mu_s),2);
    end
    [~,idx] = sort(distances);
    n_p = 1;
    mu_mid_scaled = [];
    mu_mid = [];
    if (m == 1) || (m == num_snaps)
        while (n_p <= (num_params)) && (n_p <= max(idx))
            mu_close = S_scaled(idx(n_p+1),:);
            mu_mid_scaled(n_p, :) = (mu_close + mu_s)/2;
            n_p = n_p + 1;
        end
    else
        while (n_p <= (num_params+1)) && (n_p <= max(idx))
            mu_close = S_scaled(idx(n_p+1),:);
            mu_mid_scaled(n_p, :) = (mu_close + mu_s)/2;
            n_p = n_p + 1;
        end
    end
    for i=1:num_params
        mu_min = min(S_params(:,i));
        mu_max = max(S_params(:,i));
        for j = 1:length(mu_mid_scaled)
            mu_mid(j, i) = mu_mid_scaled(j, i)*(mu_max-mu_min) + mu_min;
        end
    end
    for s = 1:length(mu_mid)
        if isempty(P_mu)
            P_mu = [P_mu, mu_mid(s)];
        elseif num_params == 1
            if ~(ismember(mu_mid(s), P_mu)) && (~(ismember(mu_mid(s), S_params)))
               P_mu = [P_mu, mu_mid(s)];
            end
        else
            if (~(ismember(mu_mid(s), P_mu, "rows"))) && (~(ismember(mu_mid(s), S_params, "rows")))
                P_mu = [P_mu, mu_mid(s)];
            end
        end
    end
end

P_mu = P_mu.';
        
end

