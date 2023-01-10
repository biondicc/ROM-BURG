function [V_H, w_ref, max_error, S_mu, P_mu] = uniform_sampling(tol,S_mu,N,x,ic,num_snaps)
% UNIFORM_SAMPLING build a basis via uniform parameter sampling with a
% specified number of snapshots
%
% INPUTS:
% <tol> = error tolerance
% <S_mu> = matrix of starting snapshot parameter locations
% <N> = dimension of the FOM
% <x> = vector of x-values
% <ic> = initial condition for the PDE
% <num_snaps> = number of snapshots to use in ROM
%
% OUTPUTS:
% <V_H> = final basis for the ROM
% <w_ref> = reference solution (average of snapshots)
% <max_error> = last maximum error value from the RBF
% <S_mu> = final matrix of snapshot parameter locations
% <P_mu> = parameter locations for error estimates/building RBF

S_size = size(S_mu);
num_params = S_size(2);
S_new = zeros(num_snaps, num_params);
for i=1:num_params
    mu_min = min(S_mu(:,i));
    mu_max = max(S_mu(:,i));
    S_new(:,i) = linspace(mu_min,mu_max,num_snaps).';
end
S_mu = S_new;

S_w = zeros(N, length(S_mu));
for i=1:length(S_mu)
    S_w(:,i) = steady_state_solver(N,S_mu(i));
end

[V_H, w_ref] = compute_POD_basis(S_w);
[P_mu] = ROM_params(S_mu);
%P_mu = P_mu_adapt;
ep_mu = [];
for j=1:length(P_mu)
    [w_bar, ~] = LSPG(V_H, w_ref, tol, x, P_mu(j), ic);
    [J_error] = error_ROM_FOM(w_bar, x, P_mu(j), ic);
    ep_mu = [ep_mu, J_error];
end
[~, max_error] = max_error_param(S_mu, P_mu, ep_mu.');


end

