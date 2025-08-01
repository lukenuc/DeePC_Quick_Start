
% =========================================================================
% FUNCTION NAME: interiorPoint.m
% AUTHOR:        Luke Nuculaj
% DESCRIPTION:   Takes the matrices for a finite-horizon OCP problem of
% interest as input and performs several interior point iterations,
% tightening the relaxation parameter \mu with each step. Returns the
% approximate optimal control input arrived at.
%
% INPUTS:
%   - H : objective function Hessian matrix
%   - f : objective function gradient
%   - A: equality constraint gradient
%   - b: equality constraint bound
%   - A_i: inequality constraint gradient
%   - b_i: inequality constraint bound
%   - alg: (1) for GE, (2) for SOR, (3) for Schur-Newton, other for A\b
%
% OUTPUTS:
%   - u_mpc: optimal control move
%
% EXIT FLAGS:
%   - exitflag: = 0 means the function was successful
%
% =========================================================================

function [u_mpc, exitflag] = interiorPoint(H, f, A, b, A_i, b_i, n_u, N, alg)

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
s = init_guess(1:ns);
S = diag(s);
u = init_guess(ns+1:end);

% anonymous functions
c_i = @(x) A_i*x - b_i;

while mu > 0.1

    H_kkt = [H zeros(size(H,1), size(Z, 2)) A_i';
        A_i eye(size(A_i, 1), size(Z, 2)) zeros(size(A_i,1), size(A_i,1));
        zeros(size(Z,1), size(H,2)) Z S];
    grad_kkt = [H*u-(A_i'*z)+f; c_i(u) + s; S*z - mu.*ones(size(S,1),1)];
    n = length(grad_kkt); 

    E = @(u, s, z) max([norm(H*u-(A_i'*z)+f, inf) ...
        norm(diag(s)*z - mu.*ones(length(s),1), inf) ...
        norm(c_i(u) + s,inf)]);
    while (E(u, s, z) > mu)

        if alg == 1 % Gaussian elimination
            [U,c,~] = GaussElim_IP(H_kkt, -grad_kkt, N, m); 
            [p,~] = backSubs_IP(U, c, N, m);
        elseif alg == 2 % SOR
            w = 0.45;
            [p,~,~] = SOR_IP(H_kkt, -grad_kkt, zeros(n,1), w, TOL, MAXIT, N, m);
        elseif alg == 3 % Schur-Newton
            inv_H_kkt = schurInverseNewton(H_kkt, N, m, TOL, MAXIT);
            p = -inv_H_kkt*grad_kkt;
        else % MATLAB built-in
            p = -H_kkt\grad_kkt; 
        end
        
        lu = length(u); ls = length(s); lz = length(z);
        p_u = p(1:lu);
        p_s = p(lu+1:lu+ls);
        p_z = p(lu+ls+1:lu+ls+lz);
        alpha_s = LineSearch(s, p_s, tau);
        alpha_z = LineSearch(z, p_z, tau);
        u = u + 50*alpha_s*p_u;
        s = s + 50*alpha_s*p_s;
        z = z + 50*alpha_z*p_z;
        S = diag(s); Z = diag(z);

        H_kkt = [H zeros(size(H,1), size(Z, 2)) A_i';
            A_i eye(size(A_i, 1), size(Z, 2)) zeros(size(A_i,1), size(A_i,1));
            zeros(size(Z,1), size(H,2)) Z S];

        grad_kkt = [H*u-(A_i'*z)+f; c_i(u) + s; S*z - mu.*ones(size(S,1),1)];

    end
    mu = sigma*mu; % tighten mu
end
u_mpc = u(1:n_u);
end

