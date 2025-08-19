
clear all; 
close all; 
rng(0); 

%% System Dynamics

m = 1; p = 1;  
Ts = 0.1; 
b = 0.1; % drag coefficient
c = 0.01; 
A = [0 1; -c -b];
B = [0; 1];
C = [1 0];
D = 0;

% discrete time dynamics
Ad = eye(2) + Ts.*A;
Bd = Ts.*B;

%% Generate persistently exciting trajectory

tf = 50;
N = tf/Ts + 1;        % number of samples in trajectory
L = 10;         % prediction horizon (referred to as "lag" in Willems' paper)
n = 2;          % state cardinality (how many states are in the system)
walk = cumsum(rand(1,N) - 0.5);
% walk = 10*sin(0.5.*(1:N).*sin(1:N)).*sin(0.05.*(1:N)) + 0.5*rand(1,N) - 0.25; 
inp = @(t_in) walk(t_in);
dynamics = @(x,u) Ad*x + Bd*u; 
measurement = @(x, u) C*x + D*u; 
x = [0; 0];     % initial state
t = linspace(0, tf, N);
u_arr = zeros(1, N); 
x_arr = zeros(2, N); 

% get trajectories and plot
for i = 1:N
    x = dynamics(x, inp(i)); 
    x_arr(:,i) = x; 
    u_arr(:,i) = inp(i); 
end

plot(t, x_arr(2,:), 'LineWidth', 2)
hold on 
plot(t, x_arr(1,:), 'LineWidth', 2)
hold on
plot(t, u_arr, 'LineWidth', 2)
hold off
legend('Vel', 'Pos', 'Accel (u)')

% generate Hankel matrices
Hy = zeros(p*(L+n), N-L+1); Hu = zeros(m*(L+n), N-L+1); 
for i = 1:N-L-n+1
    u_t = u_arr(:,i:i+m*(L+n)-1); x_t = x_arr(:,i:i+L+n-1);
    y_t = kron(eye(L+n), C)*x_t(:); 
    Hu(:,i) = u_t(:); 
    Hy(:,i) = y_t(:);
end

fprintf("Rank of %d-row Hankel matrix: %d.", size([Hu; Hy], 1), rank([Hu; Hy]));