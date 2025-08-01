% This function performs backward substitution on an upper triangular 
% system Ux = c.
% Can call this function as [x,exitflag] = backSubs(U,c)
% It is assumed that U is in upper triangular form, so Gaussian elimination
% needs to be performed before this code.
% Inputs:
%   - U: square nxn upper triangular matrix
%   - c: nx1 rhs vector
% Outputs:
%   - x: solution vector
%   - exitflag: = -1 means U is not square, 
%               = -2 means c is not a column vector or size not 
%                    compatible with U,
%               = 0 means the algorithm was successful.
%
function [x,exitflag] = backSubs(U,c)

n = size(U,1); % number of rows
m = size(U,2); % number of columns

% Check that U is square
if n ~= m
    exitflag = -1;
    return
end

% Check c is column vector and compatible with U
if n ~= size(c,1) || size(c,2) ~= 1
    exitflag = -2;
    return
end

x = zeros(n,1);

x(n) = c(n)/U(n,n);

for i = n-1 : -1 : 1
    s = 0;
    for j = i+1 : n
        s = s + U(i,j)*x(j);
    end
    x(i) = (c(i) - s)/U(i,i);
end

exitflag = 0;