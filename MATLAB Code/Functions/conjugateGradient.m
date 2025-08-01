% =========================================================================
% FUNCTION NAME: conjugateGradient.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:  Performs conjugate gradient on the SPD system. 
%
% INPUTS:
%   - A: nxn matrix
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

function [x, k] = conjugateGradient(A, b, M, x0, TOL, MAXIT)

if (isempty(M))
    M = ones(length(b), 1); 
end
x = x0(:); 
r = A*x0 - b; 
y = M.*r; 
p = -y; 

for k = 1:MAXIT
    if norm(r, inf) < TOL
        return
    else
        alpha = r'*y/(p'*A*p); 
        x = x + alpha*p; 
        rp1 = r + alpha*A*p;
        yp1 = M.*rp1; 
        beta = rp1'*yp1/(r'*y); 
        p = -yp1 + beta*p; 
        r = rp1; y = yp1; 
    end
end

end

