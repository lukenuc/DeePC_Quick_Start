% This function performs forward substitution on an upper triangular 
% system Ly = b.
% Can call this function as [y,exitflag] = forwardSubs(L,b)
% It is assumed that L is in upper triangular form, so Gaussian elimination
% needs to be performed before this code.
% Inputs:
%   - L: square nxn upper triangular matrix
%   - b: nx1 rhs vector
% Outputs:
%   - y: solution vector
%   - exitflag: = -1 means L is not square, 
%               = -2 means b is not a column vector or size not 
%                    compatible with L,
%               = 0 means the algorithm was successful.
%
function [y,exitflag] = forwardSubs(L,b)

n = size(L,1); % number of rows
m = size(L,2); % number of columns

% Check that L is square
if n ~= m
    exitflag = -1;
    return
end

% Check b is column vector and compatible with L
if n ~= size(b,1) || size(b,2) ~= 1
    exitflag = -2;
    return
end

y = zeros(n,1);

y(1) = b(1)/L(1,1);

for i = 1 : n
    s = 0;
    for j = 1 : i-1
        s = s + L(i,j)*y(j);
    end
    y(i) = (b(i) - s)/L(i,i);
end

exitflag = 0;