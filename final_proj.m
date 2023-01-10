%% Initial Problem Setup

N = 1024;
x_dom = [0, 100];
dx = (x_dom(2)-x_dom(1))/N;
x = (1:N)*dx;

ic = 1;

%% Adaptive Sampling

S_mu = linspace(0.01, 0.1, 3);
S_mu = S_mu.';
tol = 10^-4;
[V_adapt, w_ref_adapt, e_max_adapt, S_mu_adapt, P_mu_adapt] = goal_adaptive_sampling(tol,S_mu,N,x,ic);

size_V = size(V_adapt);
num_snaps = size_V(2);
%[V_uni, w_ref_uni, e_max_uni, S_mu_uni, P_mu_uni] = uniform_sampling(tol,S_mu,N,x,ic,num_snaps);

%% Uniform Sampling

[V_uni, w_ref_uni, e_max_uni, S_mu_uni, P_mu_uni] = uniform_sampling(tol,S_mu,N,x,ic,num_snaps);

%% Plotting and Errors

%b_test = 0.015;
b_test = 0.095;
w_exact = steady_state_solver(N,b_test);
J_exact = obj_func(w_exact,dx);

w_b_adapt = LSPG(V_adapt, w_ref_adapt, tol, x, b_test, ic);
w_b_uni = LSPG(V_uni, w_ref_uni, tol, x, b_test, ic);

L_error_adapt = max(abs(w_b_adapt - w_exact));
L_error_uni = max(abs(w_b_uni - w_exact));

J_adapt = obj_func(w_b_adapt,dx);
J_uni = obj_func(w_b_uni,dx);

J_error_adapt = abs(J_exact-J_adapt);
J_error_uni = abs(J_exact-J_uni);

% initial 
figure(1)
hAxes = axes('NextPlot','add',...           %# Add subsequent plots to the axes,
             'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
             'XLim',[0.01 0.1],...               %#   set the x axis limit,
             'YLim',[0 eps],...             %#   set the y axis limit (tiny!),
             'Color','none');               %#   and don't use a background color
data_1 = plot(S_mu_adapt(1:4),0,'r*','MarkerSize',10);%# Plot data set 1
hold on;
data_2 = plot(P_mu_adapt(1:4),0,'b.','MarkerSize',10);  %# Plot data set 2
legend([data_1(1),data_2(1)],'Snapshot', 'Error Eval Point','Location','Southoutside')

% final
figure(2)
hAxes = axes('NextPlot','add',...           %# Add subsequent plots to the axes,
             'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
             'XLim',[0.01 0.1],...               %#   set the x axis limit,
             'YLim',[0 eps],...             %#   set the y axis limit (tiny!),
             'Color','none');               %#   and don't use a background color
data_1 = plot(S_mu_adapt,0,'r*','MarkerSize',10);%# Plot data set 1
hold on;
data_2 = plot(P_mu_adapt,0,'b.','MarkerSize',10);  %# Plot data set 2
legend([data_1(1),data_2(1)],'Snapshot', 'Error Eval Point','Location','Southoutside')

figure(3)
hAxes = axes('NextPlot','add',...           %# Add subsequent plots to the axes,
             'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
             'XLim',[0.01 0.1],...               %#   set the x axis limit,
             'YLim',[0 eps],...             %#   set the y axis limit (tiny!),
             'Color','none');               %#   and don't use a background color
data_1 = plot(S_mu_uni,0,'r*','MarkerSize',10);%# Plot data set 1
hold on;
data_2 = plot(P_mu_uni,0,'b.','MarkerSize',10);  %# Plot data set 2
legend([data_1(1),data_2(1)],'Snapshot', 'Error Eval Point','Location','Southoutside')

set(0,'defaulttextinterpreter','latex')
figure(4)
plot(x,w_exact)
hold on;
plot(x, w_b_adapt, '--')
plot(x, w_b_uni, '-.')
title('Solution for b = 0.095 for Two Snapshot Methods', 'FontSize', 14)
xlabel('x', 'FontSize', 12);
ylabel('\textbf{w}','FontSize', 12);
legend('Full-Order Solution','Goal-Oriented','Uniform','Location','best','Interpreter','latex','FontSize', 12)
saveas(gcf, 'Figures/fin.png');


