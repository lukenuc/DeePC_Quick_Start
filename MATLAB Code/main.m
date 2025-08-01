
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

selectSPD = 1; % (0) for primal-dual system, (1) for SPD system

%% MPC setup

x = zeros(2,1); % state vector
Ts = 0.1; % sampling period [s]
Tsim = 10; % simulation period [s]
t = 0:Ts:Tsim;

b = 0.05; % drag coefficient
A = [0 1; 0 -b];
B = [0; 1];
C = eye(2);
D = [0; 0];

% discrete time dynamics
Ad = eye(2) + Ts.*A;
Bd = Ts.*B;

nx = 2; % number of states
nu = 1; % number of control inputs
p = 20; % size of prediction horizon
x0_t = [0; 0]; % initial state

xmin = repmat([-100; -30], p, 1);
xmax = repmat([100; 22], p, 1);
umin = repmat([-30], p, 1);
umax = repmat([30], p, 1);
M_ab = zeros(nx*p, nu*p);
M_ak = zeros(nx*p, nx);
qy_weight = [50; 1];
qu_weight = 1;
Qy = diag(qy_weight);
Qu = diag(qu_weight);
My = diag(repmat(qy_weight, p, 1));
Mu = diag(repmat(qu_weight, p, 1));
Mc = kron(eye(p), C);

% assembling the dense matrices
for i = 1:p
    M_ak((i-1)*nx+1:(i-1)*nx+nx, 1:nx) = C*Ad^i;
    for j = 1:p
        if (j > i)
            M_ab((i-1)*nx+1:(i-1)*nx+nx, (j-1)*nu+1:(j-1)*nu+nu) = zeros(nx,nu);
        else
            M_ab((i-1)*nx+1:(i-1)*nx+nx, (j-1)*nu+1:(j-1)*nu+nu) = C*(Ad^(i-j))*Bd;
        end
    end
end

A_i = [M_ab; -M_ab; eye(nu*p); -eye(nu*p)];
H = (1/2).*((M_ab')*My*M_ab + Mu); 

% Reference Bezier curve for position tracking
t_bez = linspace(0, 1, Tsim/Ts);
t_bez = [t_bez repmat(t_bez(end), 1, p)]; % extend ending for prediction horizon
bezier_ref = kron((1-t_bez).^3, 0) + kron(3*(1-t_bez).^2.*t_bez, -5) + kron(3*(1-t_bez).*t_bez.^2, 105) + kron(t_bez.^3, 100); 

avg_time = 0; 

for i = 1:length(t)-1
    i
    pos_ref = bezier_ref(i:i+p-1)'; 
    y_ref = zeros(nx*p, 1); 
    y_ref(1:2:end) = pos_ref; 
    f = (M_ab')*My*(M_ak*x0_t - y_ref);
    b_i = [xmax-M_ak*x0_t+y_ref; -xmin+M_ak*x0_t-y_ref; umax; -umin];
    tic
    if selectSPD == 1
        [u,~] = interiorPointDense(H,f,[],[],A_i,b_i,nu,p,1);
    else
        [u,~] = interiorPoint(H,f,[],[],A_i,b_i,nu,p,1);
    end
    avg_time = avg_time + toc/length(t); 
    x0_t = Ad*x0_t + Bd*u; % update state with optimal control
    mpc_states(:, i+1) = x0_t; 
    u_mpc(i) = u; 
end


%% plot results 

figure('Position', [100 100 1200 400]); 
subplot(1, 3, 1)
stairs(t(1:end-1), u_mpc, 'LineWidth', 2); 
title('Vehicle Control: MPC Input'); 
xlabel('Time [s]')
ylabel('Applied Acceleration [m/s/s]')
grid on 
grid minor
subplot(1, 3, 2)
plot(t(1:end-1), mpc_states(1,1:end-1), 'LineWidth', 2, 'Color', 'g'); 
hold on
plot(t(1:end-1), bezier_ref(1:length(t)-1), 'b--')
hold off
title('Vehicle Control: Position'); 
xlabel('Time [s]')
ylabel('Position [m]')
legend('Actual', 'Reference', 'Location', 'northwest')
grid on 
grid minor
subplot(1, 3, 3)
plot(t(1:end-1), mpc_states(2,1:end-1), 'LineWidth', 2, 'Color', 'r'); 
title('Vehicle Control: Velocity'); 
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid on 
grid minor
