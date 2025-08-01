
% =========================================================================
% FUNCTION NAME: mtx_inv_newton.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:  Computes the inverse of the provided matrix using
% quasi-Newton iterations. 
%
% INPUTS:
%   - A: square nxn matrix
%   - TOL: tolerance 
%   - MAXIT: maximum number of iterations
%   - x0: initial guess
%
% OUTPUTS:
%   - x: computed inverse
%   - iterates: array of all the quasi-Newton iterates
%
% EXIT FLAGS: (none)
% 
% =========================================================================

function [x, iterates] = mtx_inv_newton(A, TOL, MAXIT, x0)

n = size(A,1); 
R = @(x, A) eye(n) - A*x; 

x = x0;
for i = 1:MAXIT
    if norm(R(x, A), inf) < TOL
        return
    else
        x = x + x*R(x, A); 
    end
end

