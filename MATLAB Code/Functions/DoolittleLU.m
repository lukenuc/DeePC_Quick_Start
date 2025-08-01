
% =========================================================================
% FUNCTION NAME: DoolittleLU.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:  Performs Doolittle LU decomposition on a matrix A, 
% overwriting the contents of A with the lower and upper triangular
% matrices.
%
% INPUTS:
%   - A: square nxn matrix
%
% OUTPUTS:
%   - L: lower triangular matrix of the LU decomposition
%   - U: upper triangular matrix of the LU decomposition
%
% EXIT FLAGS: (none)
% 
% =========================================================================

function [L,U] = DoolittleLU(A)

n = size(A, 1); 
for k = 1:n

    % (1)
    sum = 0; 
    for s = 1:k-1
        sum = sum + A(k,s)*A(s,k); 
    end
    A(k,k) = A(k,k) - sum; 
    
    % (2)
    for j = k+1:n
        sum = 0; 
        for s = 1:k-1
            sum = sum + A(k,s)*A(s,j);
        end
        A(k,j) = A(k,j) - sum; 
    end

    % (3)
    for i = k+1:n
        sum = 0; 
        for s = 1:k-1
            sum = sum + A(i,s)*A(s,k);
        end
        A(i,k) = (A(i,k) - sum)/A(k,k); 
    end
    
end
L = tril(A, -1) + eye(n); 
U = triu(A); 
end

