
% =========================================================================
% FILE NAME:    main.m
% AUTHOR:       Luke Nuculaj
%
% Description:  The following script takes in Hankel matrices, cost
% function, and system constraints, and solves for the optimal control for
% the associated DeePC problem. 
% =========================================================================

clear all; close all; clc

%% System initialization

m = 1;   % No. of inputs
n = 2;   % No. of states 
p = 1;   % No. of outputs
L = 10;   % prediction horizon
N = 50;   % trajectory length
Hu = randn((n+L)*m, N-L+1);  % input Hankel matrix
Hy = randn((n+L)*p, N-L+1);  % output Hankel matrix
Hz = [Hu; Hy]; % combined Hankel matrix 

assert(L <= N, 'Prediction horizon must be less than trajectory length');
assert( size(Hu,1)==(n+L)*m, ...
    'Input Hankel matrix must have %d rows, but has %d.', ...
    (n+L)*m, size(Hu,1) );
assert( size(Hu,2)==N-L+1, ...
    'Input Hankel matrix must have %d columns, but has %d.', ...
    N-L+1, size(Hu,2) );
assert( size(Hy,1)==(n+L)*p, ...
    'Output Hankel matrix must have %d rows, but has %d.', ...
    (n+L)*p, size(Hy,1) );
assert( size(Hy,2)==N-L+1, ...
    'Output Hankel matrix must have %d columns, but has %d.', ...
    N-L+1, size(Hy,2) );
if rank(Hz) < size(Hz, 1)
    warning(['Hankel matrix NOT persistently exciting. No. of rows = %d, ' ...
        'Row rank = %d.'], size(Hz,1), rank(Hz)); 
end
fprintf('SUCCESS: Hankel matrices have been properly initialized.\n'); 


%% DeePC Initialization

Ts = 0.1; % sampling period [s]
Tsim = 10; % simulation time [s]

qy_weight = 1e0; % Tune these weights accordingly
qu_weight = 1e0;
qa_weight = 1e0; 
Qy = diag(repmat(qy_weight, 1, p));
Qu = diag(repmat(qu_weight, 1, m));
My = kron(eye(L+n), Qy);
Mu = kron(eye(L+n), Qu);
Ma = qa_weight.*eye(N-L+1);  
Hu_n = Hu(1:n*m, :);
Hy_n = Hy(1:n*p, :);
H_n  = [Hu_n; Hy_n]; 
Hu_L = Hu(n*m+1:end, :); 
Hy_L = Hy(n*p+1:end, :); 
H_L  = [Hu_L; Hy_L]; 
Mu_L = Mu(n*m+1:(L+n)*m, n*m+1:(L+n)*m); 
My_L = My(n*p+1:(L+n)*p, n*p+1:(L+n)*p);
M_L  = blkdiag(Mu_L, My_L);

H = H_L'*M_L*H_L + Ma; % Hessian 
u_ref = repmat(0, m*L, 1); % TODO: Enter input reference
y_ref = repmat(0, p*L, 1); % TODO: Enter output reference
f = -Hu_L'*Mu_L*u_ref - Hy_L'*My_L*y_ref;

% TODO: Enter your constraints here in terms of the vector z = [u_{[0,L-1]};y_{[0,L-1]}]; 
Az_eq = []; 
bz_eq = [];
Az = []; 
bz = [];  

if ~isempty(Az_eq)
    Aeq = [H_n; Az_eq*H_L]; % Do NOT remove H_n from equality constraint!
else
    Aeq = H_n; % Do NOT remove H_n from equality constraint!
end

if ~isempty(Az)
    A = Az*H_L; 
else
    A = []; 
end 
b = bz;

%% Control Loop

ctrl_seq = zeros(m, length(0:Ts:Tsim)); 
for i = 0:Ts:Tsim
    z_n = []; % TODO: Enter past 'n' I/O samples here (i.e., zn = [u_{[-n,-1]}; y_{[-n,-1]}]).
    beq = [z_n; bz_eq]; % Do NOT remove z_n from equality constraint! 
    alpha_opt = quadprog(H, f, A, b, Aeq, beq);
    u_opt = [zeros(m, m*n) eye(m, m) zeros(m, m*(L-1))]*Hu*alpha_opt;

    % TODO: Apply u_opt to your system dynamics. 
    
    ctrl_seq(:,i) = u_opt;
end