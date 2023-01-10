function [UU_back] = steady_state_solver(n,mu)
% STEADY_STATE_SOLVER solve for the steady state solution of the Burgers'
% equation using the upwind method
%
% INPUTS:
% <n> = number of grid points
% <mu> = parameter values
%
% OUTPUTS:
% <UU_back> = solution from upwind method

% discretization parameters
x_end = 100;
 
% setup spatial grid
dx = (x_end-0)/n;
xx = (1:n)*dx;

UU = zeros(n,1);

% record the initial state
UU_back = UU;

for N = 1:n
    % Upwind
    if N == 1
        UU_back_sqr =  (1)^2 + 2*dx*source(xx(N),mu);
        UU_back(N) =  sqrt(UU_back_sqr);
    else
        UU_back_sqr =  (UU_back(N-1))^2 + 2*dx*source(xx(N),mu);
        UU_back(N) =  sqrt(UU_back_sqr);
    end
end

end

