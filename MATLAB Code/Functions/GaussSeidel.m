% This function uses GaussSeidel's method to solve a linear system with 
% matrix A and rhs b.
% Can call this function as [x,k,exitflag] = GaussSeidel(A,b,x0,TOL,MAXIT)
% Inputs:
%   - A:     square nxn matrix
%   - b:     nx1 rhs vector
%   - x0:    initial guess
%   - TOL:   tolerance
%   - MAXIT: maximum number of iterations
% Outputs:
%   - x: solution (assuming the method converged)
%   - k: number of iterations
%   - exitflag: = -1 the method did not converge,
%               =  0 the method converged.
%
function [x,k,exitflag] = GaussSeidel(A,b,x0,TOL,MAXIT)

x = x0;
n = length(x);

for k = 1 : MAXIT
    % Save previous iteration
    y = x;

    for i = 1 : n
        s = 0;
        for j = 1 : n
            if i == j
                continue;
            end
            s = s + A(i,j)*x(j); % Use current iteration whenever possible
        end
        x(i) = (b(i) - s)/A(i,i);
    end
    
    if norm(x-y,'inf') < TOL
        exitflag = 0;
        return;
    end
end

exitflag = -1;