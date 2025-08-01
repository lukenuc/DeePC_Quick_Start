% This function uses SOR's method to solve a linear system with 
% matrix A and rhs b.
% Can call this function as [x,k,exitflag] = SOR(A,b,x0,omega,TOL,MAXIT)
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
function [x,k,exitflag] = SOR(A,b,x0,omega,TOL,MAXIT)

x = x0(:);
n = length(b);

for k = 1 : MAXIT
    y = x; % Save previous iteration

    for i = 1 : n
        sum1 = 0;
        for j = 1 : i-1
            sum1 = sum1 + A(i,j)*x(j);
        end
        sum2 = 0;
        for j = i+1 : n
            sum2 = sum2 + A(i,j)*y(j);
        end        
        if i == 141 
            debug = 1; 
        end
        x(i) = (omega*b(i) - omega*sum1 - omega*sum2 + (1-omega)*A(i,i)*y(i))/A(i,i);
    end

    if norm(y-x,inf) < TOL
        exitflag = 0;
        return
    end
end

exitflag = -1;