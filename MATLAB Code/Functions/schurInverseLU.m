
% =========================================================================
% FUNCTION NAME: schurInverseLU.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:  This function computes the inverse of the primal-dual
% system matrix via the Schur complement, computing the smaller inverse via
% LU factorization. 
%
% INPUTS:
%   - M: square nxn matrix
%   - p,q: dimensions
%
% OUTPUTS:
%   - Minv: inverse of M
%
% EXIT FLAGS: (none)
% =========================================================================

function [Minv] = schurInverseLU(M, p, q)

persistent Dinv; 
persistent A; 
persistent B; 
persistent C; 

if isempty(A)
    A = M(1:p,1:p); 
    B = M(1:p,p+1:end);
    C = M(p+1:end,1:p);
    Dinv = [eye(q) zeros(q); zeros(q, 2*q)];
end

D = M(p+1:end,p+1:end);

Z = D(q+1:end, 1:q); S = D(q+1:end, q+1:end);
diag_sinv = 1 ./ diag(S);
diag_z = diag(Z);

for i = 1:q
    Dinv(i+q,i) = -diag_sinv(i)*diag_z(i); 
    Dinv(i+q,i+q) = diag_sinv(i); 
end

MsD = A - B*Dinv*C;
MsDinv = mtx_inv_LU(MsD); 
T = MsDinv*B*Dinv; U = Dinv*C;
Minv = [MsDinv -T; -U*MsDinv Dinv+U*T];

end

