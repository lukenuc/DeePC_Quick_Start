
% =========================================================================
% FILE NAME:    main.m
% AUTHOR:       Luke Nuculaj
% CLASS:        APM 5334 - Applied Numerical Methods
%
% DESCRIPTION:  This script sets up the control problem and all of the
% associated matrices, which are then fed to the interior point functions
% for solving. In the for loop, select between "interiorPoint" (primal-dual
% system) and "interiorPointDense" (SPD system). The results of the
% repeated interior point are plotted, showing the control expenditure, the
% reference tracking performance, and the speed of the vehicle. 
% =========================================================================

clear all; close all; clc

run getHankel.m

%% MPC setup

x = zeros(2,1); % state vector
Ts = 0.1; % sampling period [s]
Tsim = 10; % simulation period [s]
t = 0:Ts:Tsim;

nx = 2; % number of states
nu = 1; % number of control inputs
x0_t = [8; 0]; % initial state
y_init = zeros(p*n, 1); 
y_target = repmat([10; 0], n, 1); 
u_init = zeros(m*n, 1);
u_target = zeros(m*n, 1); 

s = 2; 
qy_weight = 2e1;
qu_weight = 2e0;
qa_weight = 1e-0; 
qs_weight = 1e4; 
Qy = diag([qy_weight qy_weight]);
Qu = diag(qu_weight);
Qa = diag(qa_weight); 
Qs = diag([10*qs_weight qs_weight]); 
My = kron(eye(L+n), Qy);
Mu = kron(eye(L+n), Qu);
Ma = kron(eye(N-L+1), Qa); 
Ms = kron(eye(L+n), Qs); 
H = blkdiag(Mu, My, Ma, Ms); 

% assembling the dense matrices
q = (m + p + s)*(L + n) + (N - L + 1);
nm = m*(L+n); np = p*(L+n); na = N - L + 1; ns = s*(L+n);  
Km_left = [eye(m*n) zeros(m*n, m*L)]; 
Kp_left = [eye(p*n) zeros(p*n, p*L)];
Km_right = [zeros(m*n, m*L) eye(m*n)];
Kp_right = [zeros(p*n, p*L) eye(p*n)];

t_bez = linspace(0, 1, Tsim/Ts);
t_bez = [t_bez repmat(t_bez(end), 1, p)]; % extend ending for prediction horizon
bezier_ref = 1 + kron((1-t_bez).^3, 8) + kron(3*(1-t_bez).^2.*t_bez, 7.5) + kron(3*(1-t_bez).*t_bez.^2, 10.5) + kron(t_bez.^3, 10); 
% 
% A_eq = [eye(nm) zeros(nm, np) -Hu zeros(nm, ns); ...
%         zeros(np,nm) eye(np) -Hy eye(ns); ...
%         Km_left zeros(m*n, q-nm);
%         zeros(p*n, nm) Kp_left zeros(p*n, na+ns);
%         Km_right zeros(m*n, q-nm);
%         zeros(p*n, nm) Kp_right zeros(p*n, na + ns)]; 
% u_mpc = zeros(m, length(t)); 
% mpc_states = zeros(p, length(t));
% mpc_states(:,1) = x0_t;
% u_init = bezier_ref(1:n);
% y_init = getTrajectory(dynamics, measurement, u_init, x0_t);

%% Sparse QP Implementation
% for i = 1:length(t)-1
%     b_eq = [zeros(nm, 1); zeros(np, 1); u_init(:); y_init(:); u_target; y_target];
%     z_target = [repmat(0, L+n, 1); repmat([10; 0], L+n, 1); zeros(na, 1); zeros(ns, 1)];
%     f = (-z_target'*H)'; 
%     options = optimoptions('quadprog','Algorithm','active-set');
%     z_opt = quadprog(H, f, [], [], A_eq, b_eq, [], [], zeros(q, 1), options);
%     u = z_opt(n+1); 
%     x0_t = Ad*x0_t + Bd*u; % update state with optimal control
%     mpc_states(:, i+1) = x0_t;
%     u_mpc(:, i) = u; 
%     u_init = [u_init(2:end) u];
%     y_init = [y_init(:,2:end) C*x0_t];
% end

%% Dense Control Law

u_init = bezier_ref(1:n);
y_init = getTrajectory(dynamics, measurement, u_init, x0_t);
u_mpc = zeros(m, length(t)); 
mpc_states = zeros(n, length(t));
mpc_states(:,1) = x0_t;

Hu_n = Hu(1:n*m, :);
Hy_n = Hy(1:n*p, :);
H_n  = [Hu_n; Hy_n]; 
Hu_L = Hu(n*m+1:end, :); 
Hy_L = Hy(n*p+1:end, :); 
H_L  = [Hu_L; Hy_L]; 
Mu_L = Mu(n*m+1:(L+n)*m, n*m+1:(L+n)*m); 
My_L = My(n*p+1:(L+n)*p, n*p+1:(L+n)*p);
M_L  = blkdiag(Mu_L, My_L); 

H = H_L'*M_L*H_L + Ma; 
Hinv = inv(H); 
u_L = repmat(0, L, 1); y_L = repmat(10, L, 1);
f = -Hu_L'*Mu_L*u_L - Hy_L'*My_L*y_L; 
alpha_unc = -H\f; 
g = alpha_unc - Hinv*H_n'*inv(H_n*Hinv*H_n')*H_n*alpha_unc; 
gu = [zeros(m, m*n) eye(m, m) zeros(m, m*(L-1))]*Hu*g;
Pu = [zeros(m, m*n) eye(m, m) zeros(m, m*(L-1))]*Hu*Hinv*H_n'*inv(H_n*Hinv*H_n');

for i = 1:length(t)-1 
    z = [u_init(:); y_init(:)]; 
    u = gu + Pu*z; 
    x0_t = Ad*x0_t + Bd*u; % update state with optimal control
    mpc_states(:, i+1) = x0_t;
    u_mpc(:, i) = u; 
    u_init = [u_init(2:end) u];
    y_init = [y_init(2:end) C*x0_t];
end

%% plot results 

figure('Position', [100 100 1200 400]); 
subplot(1, 3, 1)
stairs(t(1:length(u_mpc)), u_mpc, 'LineWidth', 2); 
title('Vehicle Control: MPC Input'); 
xlabel('Time [s]')
ylabel('Applied Acceleration [m/s/s]')
grid on 
grid minor
subplot(1, 3, 2)
plot(t, mpc_states(1,:), 'LineWidth', 2, 'Color', 'g'); 
title('Vehicle Control: Position'); 
xlabel('Time [s]')
ylabel('Position [m]')
grid on 
grid minor
subplot(1, 3, 3)
plot(t, mpc_states(2,:), 'LineWidth', 2, 'Color', 'r'); 
title('Vehicle Control: Velocity'); 
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid on 
grid minor
