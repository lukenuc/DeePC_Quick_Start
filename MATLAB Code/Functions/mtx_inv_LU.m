
% =========================================================================
% FUNCTION NAME: mtx_inv_LU.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:  Computes the inverse of the provided matrix using
% the LU factorization and forward + backward substitution on
% consecutive linear systems constructed with the standard basis
% vectors. 
%
% INPUTS:
%   - A: square nxn matrix
%
% OUTPUTS:
%   - x: computed inverse
%
% EXIT FLAGS: (none)
% 
% =========================================================================

function [x] = mtx_inv_LU(A)

n = size(A, 1); 
[L,U] = DoolittleLU(A); 
x = eye(n); 

for i = 1:n
    ei = x(:,i); 
    [ci,~] = forwardSubs(L,ei); 
    [x(:,i),~] = backSubs(U,ci); 
end

end

