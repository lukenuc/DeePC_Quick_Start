% =========================================================================
% FUNCTION NAME: backSubs_IP.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:  Performs backward substitution optimized for the specific
% sparse structure of the interior point problem. 
%
% INPUTS:
%   - U: square nxn matrix
%   - c: vector of n elements
%   - p, q: dimension variables
%
% OUTPUTS:
%   - x: solution to the linear system
%
% EXIT FLAGS: = -2, incompatible U and C
%             = -1, U not square
%             = 0, successful
% 
% =========================================================================
function [x,exitflag] = backSubs_IP(U,c,p,q)

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

for i = n: -1 : n-q+1
    x(i) = c(i)/U(i,i);
end

for i = n-q : -1 : p+1
    s = 0;
    for j = n-q+1 : n
        s = s + U(i,j)*x(j);
    end
    x(i) = c(i) - s;
end

for i = p : -1 : 1
    s = 0;
    for j = i+1 : p
        s = s + U(i,j)*x(j);
    end
    for j = n-q+1:n
        s = s + U(i,j)*x(j);
    end
    x(i) = (c(i) - s)/U(i,i);
end

exitflag = 0;