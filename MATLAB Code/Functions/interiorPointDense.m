
% =========================================================================
% FUNCTION NAME: interiorPointDense.m
% AUTHOR:        Luke Nuculaj
% DESCRIPTION:   Takes the matrices for a finite-horizon OCP problem of
% interest as input and performs several interior point iterations,
% tightening the relaxation parameter \mu with each step. Returns the
% approximate optimal control input arrived at. Constructs the KKT linear
% system using the double augemented SPD form. User may choose between
% conjugate gradient (default) and SOR. 
%
% INPUTS:
%   - H : objective function Hessian matrix
%   - f : objective function gradient
%   - A: equality constraint gradient
%   - b: equality constraint bound
%   - A_i: inequality constraint gradient
%   - b_i: inequality constraint bound
%   - alg: (1) SOR, other for conjugate gradient
%
% OUTPUTS:
%   - u_mpc: optimal control move
%
% EXIT FLAGS:
%   - exitflag: = 0 means the function was successful
%
% =========================================================================

function [u_mpc, exitflag] = interiorPointDense(H, f, A, b, A_i, b_i, n_u, N, alg)

exitflag = 0;

% settings
mu = 100;
sigma = 0.5;
tau = 0.995;
TOL = 1.e-9; 
MAXIT = 1000; 

% initial guesses for Lagrange multipliers (y, z) and slacks (s)
y = ones(size(A, 1), 1);
z = ones(size(A_i, 1), 1);
Z = diag(z);
m = length(b_i);

% Phase 1 method for generating feasible initial guess
nx = size(H,1); ns = length(b_i);
f_init = [ones(ns, 1); zeros(nx, 1)];
A_init = [-eye(ns) zeros(ns,nx);
    eye(ns) A_i];
b_init = [-5.*ones(ns,1); b_i];
options = optimoptions('linprog', 'Display', 'off');
init_guess = linprog(f_init,A_init,b_init,[],[],[],[],options);
u = init_guess(ns+1:end);

% anonymous functions
c_i = @(x) A_i*x - b_i;

min_eigs = []; 

while mu > 0.1

    Zinv = diag(1 ./ z);
    H_kkt = [H + 2*A_i'*(1/mu)*Z*Z*A_i -A_i'; -A_i mu*Zinv*Zinv];
    grad_kkt = [H*u-(A_i'*z)+f; -c_i(u) - mu*Zinv*ones(m,1)];
    r1 = grad_kkt(1:N); r2 = grad_kkt(N+1:end);
    grad_kkt = [r1 - 2*A_i'*(1/mu)*Z*Z*r2; r2];
    n = length(grad_kkt); 
    min_eigs = [min_eigs min(eig(H_kkt))]; 

    c_i = @(x) A_i*x - b_i;
    E = @(u, z, Zinv) max([norm(H*u-(A_i'*z)+f, inf) ...
        norm(c_i(u) + mu*Zinv*ones(m,1), inf)]);
    while (E(u, z, Zinv) > mu)

        % convert to COO
        [row, col, v] = find(H_kkt);
        coo = sortrows([row col v], 1);
        row = coo(:,1); col = coo(:,2); v = coo(:,3);
        
        if alg == 1 
            % SOR
            [p,~,~] = SOR_COO(row, col, v, -grad_kkt, zeros(n,1), 1, TOL, MAXIT);
        else
            % conjugate gradient method
            [p,~] = CG_COO(row, col, v, -grad_kkt, [], zeros(n, 1), TOL, MAXIT);
        end
        
        % line search and step
        lu = length(u); lz = length(z);
        p_u = p(1:lu);
        p_z = p(lu+1:lu+lz);
        alpha_z = LineSearch(z, p_z, tau);
        u = u + 50*alpha_z*p_u;
        z = z + 50*alpha_z*p_z;
        Z = diag(z);
        Zinv = diag(1 ./ z);

        % update linear system
        H_kkt = [H + 2*A_i'*(1/mu)*Z*Z*A_i -A_i'; -A_i mu*Zinv*Zinv];
        grad_kkt = [H*u-(A_i'*z)+f; -c_i(u) - mu*Zinv*ones(m,1)];
        r1 = grad_kkt(1:N); r2 = grad_kkt(N+1:end);
        grad_kkt = [r1 - 2*A_i'*(1/mu)*Z*Z*r2; r2];

    end
    mu = sigma*mu; % tighten mu
end
u_mpc = u(1:n_u);
end

