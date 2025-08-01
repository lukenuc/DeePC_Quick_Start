
% =========================================================================
% FUNCTION NAME: CG_COO.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:  Performs conjugate gradient where the sparse matrix is
% represented in COO format. 
%
% INPUTS:
%   - arow: row index vector
%   - acol: column index vector
%   - aval: value vector
%   - b: RHS vector
%   - M: preconditioner (as a vector)
%   - x0: initial guess
%   - TOL: tolerance for convergence
%   - MAXIT: maximum number of iterations
%
% OUTPUTS:
%   - x: solution to the linear system
%   - k: iterations to converge
%
% EXIT FLAGS: = -2, incompatible U and C
%             = -1, U not square
%             = 0, successful
% 
% =========================================================================

function [x,k] = CG_COO(arow, acol, aval, b, M, x0, TOL, MAXIT)

if (isempty(M))
    M = ones(length(b), 1); 
end
x = x0(:); 
r = mvmCOO(arow, acol, aval, x) - b; 
y = M.*r; 
p = -y; 

for k = 1:MAXIT
    if norm(r, inf) < TOL
        return
    else
        t = mvmCOO(arow, acol, aval, p); 
        alpha = r'*y/(p'*t); 
        x = x + alpha*p; 
        rp1 = r + alpha*t;
        yp1 = M.*rp1; 
        beta = rp1'*yp1/(r'*y); 
        p = -yp1 + beta*p; 
        r = rp1; y = yp1; 
    end
end

end

