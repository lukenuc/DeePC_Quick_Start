
% =========================================================================
% FUNCTION NAME: mvmCOO.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:  This function performs matrix vector multiplication where
% the matrix is represented in COO format. 
%
% INPUTS:
%   - arow: row index vector
%   - acol: column index vector
%   - aval: value vector
%   - v: vector in the matrix vector multiplication
%
% OUTPUTS:
%   - y: result
%
% EXIT FLAGS: (none)
% =========================================================================

function [y] = mvmCOO(arow, acol, aval, v)

n = length(v); 
l = length(aval); 
y = zeros(n, 1); 

c = 1;
for i = 1:n
    sum = 0; 
    while (c <= l & arow(c)==i)
        col = acol(c); a = aval(c); 
        sum = sum + a*v(col); 
        c = c+1; 
    end
    y(i) = sum; 
end
end