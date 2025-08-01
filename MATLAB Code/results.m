
% =========================================================================
% FILE NAME:    results.m
% AUTHOR:       Luke Nuculaj
% CLASS:        APM 5334 - Applied Numerical Methods
%
% DESCRIPTION:  This script produces all of the results that appear in the
% report. 
% =========================================================================

clear all; clc; close all;
load('test.mat')
TOL = 1.e-9; MAXIT = 1000;

%% Testing execution time on problem of a fixed size (N = 20)

num_methods = 4;
num_trials = 1000;
times = zeros(num_methods, 1);
time_arr = zeros(num_methods, num_trials);
iters = zeros(num_methods, 1);
[row, col, v] = find(H_kkt);
coo = sortrows([row col v], 1);
row = coo(:,1); col = coo(:,2); v = coo(:,3);
p = zeros(length(grad_kkt),1);
n = length(grad_kkt);

for j = 1:num_trials
    j
    % Gaussian elimination
    tic
    [U,c,~] = GaussElim_IP(H_kkt,-grad_kkt,N,m);
    [p,~] = backSubs_IP(U,c,N,m);
    time_arr(1, j) = toc;

    % SOR
    tic
    [p,~,~] = SOR_COO(row, col, v, -grad_kkt, zeros(n,1), 0.45, TOL, MAXIT);
    time_arr(2, j) = toc;

    % Schur-Newton
    tic
    inv_H_kkt = schurInverseNewton(H_kkt, N, m, TOL, MAXIT);
    p = -inv_H_kkt*grad_kkt;
    time_arr(3, j) = toc;

    % MATLAB built-in inverse
    tic
    p = -H_kkt\grad_kkt;
    time_arr(4, j) = toc;
end

for i = 1:num_methods
    t = time_arr(i,:);
    t(t == max(t) | t == min(t)) = [];
    times(i) = mean(t);
end

% Check iterations to converge for SOR as a function of omega
w_arr = 0:0.01:2;
iters = zeros(1, length(w_arr));
for i = 1:length(w_arr)
    [~,k,~] = SOR_COO(row, col, v, -grad_kkt, zeros(n,1), w_arr(i), TOL, MAXIT);
    iters(i) = k;
end

% Plot execution times as a bar graph
figure;
subplot(1, 2, 1)
methods = {'GE', 'SOR', 'Schur-Newton', 'A\b'};
bar(times, 'b');
set(gca, 'XTickLabel', methods, 'XTick', 1:num_methods);
ylabel('Execution Time (seconds)');
ylim([0 1.2*max(times)])
subplot(1, 2, 2)
plot(w_arr, iters, 'LineWidth', 2)
grid minor
xlabel('\omega'); ylabel('Iterations')


%% Testing execution times of optimized versus unoptimized methods

num_methods = 3;
num_trials = 1000;
times_opt = zeros(num_methods, 1);
times = zeros(num_methods, 1);
time_arr_opt = zeros(num_methods, num_trials);
time_arr = zeros(num_methods, num_trials);
[row, col, v] = find(H_kkt);
coo = sortrows([row col v], 1);
row = coo(:,1); col = coo(:,2); v = coo(:,3);
n = length(grad_kkt);

for j = 1:num_trials

    % Gaussian elimination
    [U,c,~] = GaussElimPivoting(H_kkt,-grad_kkt); % unoptimized
    [p,~] = backSubs(U,c);
    time_arr(1, j) = toc;
    tic
    [U,c,~] = GaussElim_IP(H_kkt,-grad_kkt,N,m);
    [p,~] = backSubs_IP(U,c,N,m);
    time_arr_opt(1, j) = toc;

    % SOR
    tic
    [p,~,~] = SOR(H_kkt, -grad_kkt, zeros(length(grad_kkt),1), 0.45, TOL, MAXIT);
    time_arr(2, j) = toc;
    tic
    [p,~,~] = SOR_COO(row, col, v, -grad_kkt, zeros(n,1), 0.45, TOL, MAXIT);
    time_arr_opt(2, j) = toc;

    % Schur-Newton
    tic
    inv_H_kkt = schurInverseLU(H_kkt,N,m);
    p = -inv_H_kkt*grad_kkt;
    time_arr(3, j) = toc;
    tic
    inv_H_kkt = schurInverseNewton(H_kkt,N,m,TOL,MAXIT);
    p = -inv_H_kkt*grad_kkt;
    time_arr_opt(3, j) = toc;

end

for i = 1:num_methods
    t = time_arr(i,:);
    t(t == max(t) | t == min(t)) = [];
    times(i) = mean(t);

    t = time_arr_opt(i,:);
    t(t == max(t) | t == min(t)) = [];
    times_opt(i) = mean(t);
end

figure;
methods = {'GE', 'SOR', 'Schur-Newton'};
bar_data = [times, times_opt];
bar_colors = [1 0 1; 0 1 1]; % Magenta and cyan colors

% Create the bar plot
b = bar(bar_data, 'grouped');

% Set bar colors
for k = 1:numel(b)
    b(k).FaceColor = 'flat';
    b(k).CData = repmat(bar_colors(k,:), numel(methods), 1);
end

% Customize the graph
set(gca, 'XTickLabel', methods, 'XTick', 1:num_methods);
xlabel('Method');
ylabel('Execution Time (seconds)');
legend({'Unoptimized', 'Optimized'}, 'Location', 'northeast');
grid on;
grid minor;

%% Testing the symmetrized double-augmented IP method

Sigma = mu.*(inv(S))^2;
H_sym = [H + 2*A_i'*Sigma*A_i -A_i'; -A_i inv(Sigma)];
grad_kkt = [H*u-(A_i'*z)+f; c_i(u) + mu*inv(Z)*ones(m,1)];
r1 = grad_kkt(1:N); r2 = grad_kkt(N+1:end);
grad_kkt = [r1 - 2*A_i'*Sigma*r2; r2];

num_methods = 3;
num_trials = 1000;
times = zeros(num_methods, 1);
time_arr = zeros(num_methods, num_trials);
iters = zeros(num_methods, 1);
[row, col, v] = find(H_sym);
coo = sortrows([row col v], 1);
row = coo(:,1); col = coo(:,2); v = coo(:,3);
p = zeros(length(grad_kkt),1);
n = length(grad_kkt);

for j = 1:num_trials

    % SOR
    tic
    [p,~,~] = SOR_COO(row, col, v, -grad_kkt, zeros(n,1), 1, TOL, MAXIT);
    time_arr(1, j) = toc;

    % conjugate gradient method
    tic
    [p,k] = CG_COO(row, col, v, -grad_kkt, [], zeros(n, 1), TOL, MAXIT);
    time_arr(2, j) = toc;

    % MATLAB built-in inverse
    tic
    p = -H_sym\grad_kkt;
    time_arr(3, j) = toc;
end

for i = 1:num_methods
    t = time_arr(i,:);
    t(t == max(t) | t == min(t)) = [];
    times(i) = mean(t);
end

% Check iterations to converge for SOR as a function of omega
w_arr = 0:0.01:2;
iters = zeros(1, length(w_arr));
for i = 1:length(w_arr)
    [~,k,~] = SOR_COO(row, col, v, -grad_kkt, zeros(n,1), w_arr(i), TOL, MAXIT);
    iters(i) = k;
end

% Plot execution times as a bar graph
figure;
subplot(1, 2, 1)
methods = {'SOR (\omega = 1)', 'CG', 'A\b'};
bar(times, 'b');
set(gca, 'XTickLabel', methods, 'XTick', 1:num_methods);
ylabel('Execution Time (seconds)');
ylim([0 1.2*max(times)])
subplot(1, 2, 2)
plot(w_arr, iters, 'LineWidth', 2)
grid minor
xlabel('\omega'); ylabel('Iterations')

%% Optimized versus unoptimized for SPD system

Sigma = mu.*(inv(S))^2;
H_sym = [H + 2*A_i'*Sigma*A_i -A_i'; -A_i inv(Sigma)];
grad_kkt = [H*u-(A_i'*z)+f; c_i(u) + mu*inv(Z)*ones(m,1)];
r1 = grad_kkt(1:N); r2 = grad_kkt(N+1:end);
grad_kkt = [r1 - 2*A_i'*Sigma*r2; r2];
num_trials = 1000;
num_methods = 2;
times = zeros(num_methods, 1);
times_opt = zeros(num_methods, 1);
time_arr = zeros(num_methods, num_trials);
iters = zeros(num_methods, 1);
[row, col, v] = find(H_sym);
coo = sortrows([row col v], 1);
row = coo(:,1); col = coo(:,2); v = coo(:,3);
p = zeros(length(grad_kkt),1);
n = length(grad_kkt);

for j = 1:num_trials

    % SOR
    tic
    [p,~,~] = SOR(H_sym, -grad_kkt, zeros(n,1), 1, TOL, MAXIT);
    time_arr(1, j) = toc;
    tic
    [p,~,~] = SOR_COO(row, col, v, -grad_kkt, zeros(n,1), 1, TOL, MAXIT);
    time_arr_opt(1, j) = toc;

    % conjugate gradient method
    tic
    [p,k] = conjugateGradient(H_sym, -grad_kkt, [], zeros(n, 1), TOL, MAXIT);
    time_arr(2, j) = toc;
    tic
    [p,k] = CG_COO(row, col, v, -grad_kkt, [], zeros(n, 1), TOL, MAXIT);
    time_arr_opt(2, j) = toc;

end

for i = 1:num_methods
    t = time_arr(i,:);
    t(t == max(t) | t == min(t)) = [];
    times(i) = mean(t);

    t = time_arr_opt(i,:);
    t(t == max(t) | t == min(t)) = [];
    times_opt(i) = mean(t);
end

figure;
methods = {'SOR', 'CG'};
bar_data = [times, times_opt];
bar_colors = [1 0 1; 0 1 1]; % Magenta and cyan colors

% Create the bar plot
b = bar(bar_data, 'grouped');

% Set bar colors
for k = 1:numel(b)
    b(k).FaceColor = 'flat';
    b(k).CData = repmat(bar_colors(k,:), numel(methods), 1);
end

% Customize the graph
set(gca, 'XTickLabel', methods, 'XTick', 1:num_methods);
xlabel('Method');
ylabel('Execution Time (seconds)');
legend({'Unoptimized', 'Optimized'}, 'Location', 'northeast');
grid on;
grid minor;

%% See how the methods scale with prediction horizon

num_methods = 2;
num_trials = 2000;
N_arr = 5:1:30;
times = zeros(num_methods, 1);
time_arr = zeros(num_methods, length(N_arr));
iters = zeros(num_methods, 1);

for i = 1:length(N_arr)
    N = N_arr(i);
    H_sym = getFauxMatrix(N);
    [row, col, v] = find(H_sym);
    coo = sortrows([row col v], 1);
    row = coo(:,1); col = coo(:,2); v = coo(:,3);
    n = size(H_sym, 1); 
    for j = 1:num_trials

        % SOR
        tic
        [p,~,~] = SOR_COO(row, col, v, ones(n,1), zeros(n,1), 1, TOL, 1);
        time_arr(1, i) = time_arr(1, i) + toc;

        % conjugate gradient method
        tic
        [p,k] = CG_COO(row, col, v, ones(n,1), [], zeros(n, 1), TOL, 1);
        time_arr(2, i) = time_arr(2, i) + toc;

    end

    time_arr(1, i) = time_arr(1, i)/num_trials; 
    time_arr(2, i) = time_arr(2, i)/num_trials; 
end

plot(N_arr, time_arr, 'LineWidth', 2)
legend('SOR', 'CG', 'Location', 'northwest')
ylabel('Avg. Time Per Iteration (seconds)')
xlabel('Prediction Horizon')