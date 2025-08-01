
% =========================================================================
% FUNCTION NAME: GaussElim_IP.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:  This function performs Gaussian elimination on the
% primal-dual system accounting for sparsity.
%
% INPUTS:
%   - A: square nxn matrix
%   - b: nx1 rhs vector
%   - p,q: dimensions
%
% OUTPUTS:
%   - A: upper triangular form of A
%   - b: modified rhs
%
% EXIT FLAGS:
%   - exitflag  = -1 means A is not square,
%               = -2 means b is not a column vector or size not
%                    compatible with A,
%               = 0 means the algorithm was successful.
% =========================================================================

function [A,b,exitflag] = GaussElim_IP(A,b,p,q)

n = size(A,1); % number of rows
m = size(A,2); % number of columns

% Check that A is square
if n ~= m
    exitflag = -1;
    return
end

% Check b is column vector and compatible with A
if n ~= size(b,1) || size(b,2) ~= 1
    exitflag = -2;
    return
end

for j = 1 : p+q
    if (j <= p)
        for i = j+1:p+q
            m = A(i,j)/A(j,j);
            for k = j : n
                A(i,k) = A(i,k) - m*A(j,k);
            end
            b(i) = b(i) - m*b(j);
        end
    else
        m = A(j+q,j); 
        A(j+q,j) = 0;
        b(j+q) = b(j+q) - m; 
    end
end
exitflag = 0;