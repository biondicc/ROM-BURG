%% Initial Problem Setup

N = 1024;
x_dom = [0, 100];
dx = (x_dom(2)-x_dom(1))/N;
x = (1:N)*dx;

ic = 1;

%% Visualize Parameter Space

S_mu = linspace(0.01, 0.1, 3);
S_mu = S_mu.';
tol = 10^-4;

S_w = zeros(N, length(S_mu));
for i=1:length(S_mu)
    S_w(:,i) = steady_state_solver(N,S_mu(i));
end

set(0,'defaulttextinterpreter','latex')
figure(1), clf,
plot(x, S_w)
title('Snapshots of Solutions at Various Design Parameters', 'FontSize', 14)
xlabel('x', 'FontSize', 12);
ylabel('\textbf{w}','FontSize', 12);
legend('b = 0.01', 'b = 0.055', 'b = 0.1','Location','best','Interpreter','latex','FontSize', 12)
saveas(gcf, 'Figures/snaps_b_3.png');

%% Study SVD

S_mu = linspace(0.01, 0.1, 20);
S_mu = S_mu.';
tol = 10^-4;

S_w = zeros(N, length(S_mu));
for i=1:length(S_mu)
    S_w(:,i) = steady_state_solver(N,S_mu(i));
end

w_ref = mean(S_w, 2);
Q = S_w - w_ref;
[V,D,~] = svd(Q,"econ");

b_test = 0.055;
w_exact = steady_state_solver(N,b_test);
[w_bar, iter] = LSPG(V, w_ref, tol, x, b_test, ic);

sigma = D(sub2ind(size(D),1:size(D,1),1:size(D,2)));
tot_sigma = sum(sigma);
cum_sigma = zeros(1,length(sigma));
percent_sigma = zeros(1,length(sigma));
for j=1:length(sigma)
    cum_sigma(j) = sum(sigma(1:j));
    percent_sigma(j) = cum_sigma(j)/tot_sigma;
end

w_lim = zeros(N, 4);
p = 1;
for j=[1,3,4,6]
    [w_lim(:, p),~] = LSPG(V(:,1:j), w_ref, tol, x, b_test, ic);
    p=p+1;
end

plot(x,w_exact)
hold on;
plot(x, w_lim, '--')
title('Solution for b = 0.055 for Various Subspace Sizes', 'FontSize', 14)
xlabel('x', 'FontSize', 12);
ylabel('\textbf{w}','FontSize', 12);
legend('Full-Order Solution','1 Mode','3 Modes','4 Modes','6 Modes','Location','best','Interpreter','latex','FontSize', 12)
saveas(gcf, 'Figures/sub_size.png');

for m=1:4
    L_error(m) = max(abs(w_lim(:,m) - w_exact));
end

%% Two Design Parameters

mu_b = linspace(0.01, 0.1, 10);
mu_a  = linspace(0.01, 100, 10);
S_mu = zeros(length(mu_a)*length(mu_b),2);
shift = 0;
for i = 1:length(mu_b)
    b = mu_b(i);
    
    for j = 1:length(mu_a)
        S_mu(j+shift,:) = [b; mu_a(j)];
    end
    shift = shift + length(mu_a);
end

S_w = zeros(N, length(S_mu));
for i=1:length(S_mu)
    S_w(:,i) = steady_state_solver(N,S_mu(i,:));
end

w_ref = mean(S_w, 2);
Q = S_w - w_ref;
[V,D,~] = svd(Q,"econ");

sigma = D(sub2ind(size(D),1:size(D,1),1:size(D,2)));
tot_sigma = sum(sigma);
cum_sigma = zeros(1,length(sigma));
percent_sigma = zeros(1,length(sigma));
for j=1:length(sigma)
    cum_sigma(j) = sum(sigma(1:j));
    percent_sigma(j) = cum_sigma(j)/tot_sigma;
end

w_lim = zeros(N, 6);
for j=1:6
    [w_lim(:, j),~] = LSPG(V(:,1:j), w_ref, tol, x, b_test, ic);
end