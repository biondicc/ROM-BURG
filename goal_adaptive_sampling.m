function [V_H, w_ref, max_error, S_mu, P_mu] = goal_adaptive_sampling(tol,S_mu,N,x,ic)
% GOAL_ADAPTIVE_SAMPLING build a basis to a specific error tolerance by
% adding snapshots at locations with high errors
%
% INPUTS:
% <tol> = error tolerance
% <S_mu> = matrix of starting snapshot parameter locations
% <N> = dimension of the FOM
% <x> = vector of x-values
% <ic> = initial condition for the PDE
%
% OUTPUTS:
% <V_H> = final basis for the ROM
% <w_ref> = reference solution (average of snapshots)
% <max_error> = last maximum error value from the RBF
% <S_mu> = final matrix of snapshot parameter locations
% <P_mu> = parameter locations for error estimates/building RBF

S_w = zeros(N, length(S_mu));
for i=1:length(S_mu)
    S_w(:,i) = steady_state_solver(N,S_mu(i));
end

[V_H, w_ref] = compute_POD_basis(S_w);
[P_mu] = ROM_params(S_mu);

ep_mu = [];
for j=1:length(P_mu)
    [w_bar, ~] = LSPG(V_H, w_ref, tol, x, P_mu(j), ic);
    [J_error] = error_ROM_FOM(w_bar, x, P_mu(j), ic);
    ep_mu = [ep_mu, J_error];
end
[mu_max_error, max_error] = max_error_param(S_mu, P_mu, ep_mu.');

while max_error > tol
    S_w_new = steady_state_solver(N,mu_max_error);
    S_w = horzcat(S_w,S_w_new);
    [V_h, w_ref_h] = compute_POD_basis(S_w);
    ep_mu = [];
    for j=1:length(P_mu)
        [w_bar, ~] = LSPG(V_h, w_ref_h, tol, x, P_mu(j), ic);
        [J_error] = error_ROM_FOM(w_bar, x, P_mu(j), ic);
        ep_mu = [ep_mu, J_error];
    end
%     for k=1:length(P_mu)
%         [w_bar, ~] = LSPG(V_H, w_ref, tol, x, P_mu(k), ic);
%         [w_bar_h, ~] = LSPG(V_h, w_ref_h, tol, x, P_mu(k), ic);
%         [J_error] = error_coarse_fine_ROM(w_bar, w_bar_h, x, P_mu(k), ic, V_h);
%         ep_mu(k) = ep_mu(k) - J_error;
%     end
    % SKIPPING ERROR UPDATE
    P_mu_star = ROM_params_max_error(sort(S_mu),mu_max_error);
    V_H = V_h;
    w_ref = w_ref_h;
    S_mu = [S_mu; mu_max_error];
    for p=1:length(P_mu_star)
        mu_star = P_mu_star(p);
        P_mu = [P_mu; mu_star];
        [w_bar, ~] = LSPG(V_H, w_ref, tol, x, mu_star, ic);
        [J_error] = error_ROM_FOM(w_bar, x, mu_star, ic);
        ep_mu = [ep_mu, J_error];
    end
    [mu_max_error, max_error] = max_error_param(S_mu, P_mu, ep_mu.');
    if any(ismember(round(P_mu,3),round(mu_max_error,3)))
        mu_max_error = mu_max_error + 0.0005;
    end
%     for t=1:length(ep_mu)
%         if abs(ep_mu(t)) < 10^(-5)
%             return;
%         end
%     end
end

end


