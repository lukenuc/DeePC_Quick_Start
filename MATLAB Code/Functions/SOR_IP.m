
% =========================================================================
% FUNCTION NAME: SOR_IP.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:  Performs SOR on the primal-dual system where the
% iterations are optimized for the sparsity of the system matrix. 
%
% INPUTS:
%   - A: nxn matrix
%   - b: RHS vector
%   - x0: initial guess
%   - omega: relaxation parameter
%   - TOL: tolerance for convergence
%   - MAXIT: maximum number of iterations
%   - p: prediction horizon
%   - m: number of constraints
%
% OUTPUTS:
%   - x: solution to the linear system
%   - k: iterations to converge
%
% EXIT FLAGS: = -1, did not converge
%             = 0, did converge
% 
% =========================================================================

function [x,k,exitflag] = SOR_IP(A,b,x0,omega,TOL,MAXIT,p,m)

x = x0(:);
n = length(b);
t = [1:p, p+m+1:p+2*m]; 

for k = 1 : MAXIT
    y = x; % Save previous iteration

    for i = 1 : n
        if i <= p % Phase 1
            sum1 = 0;
            for j = 1 : i - 1
                sum1 = sum1 + A(i,j)*x(j);
            end
            sum2 = 0;
            for j = t(i+1:end)
                sum2 = sum2 + A(i,j)*y(j);
            end
            x(i) = (omega*b(i) - omega*sum1 - omega*sum2 + (1-omega)*A(i,i)*y(i))/A(i,i); 
        elseif i <= p+m % Phase 2
            sum1 = 0;
            for j = 1 : p
                sum1 = sum1 + A(i,j)*x(j);
            end
            x(i) = (omega*b(i) - omega*sum1 + (1-omega)*A(i,i)*y(i))/A(i,i); 
        else % Phase 3
            x(i) = (omega*b(i) - omega*A(i,i-m)*x(i-m) + (1-omega)*A(i,i)*y(i))/A(i,i); 
        end
    end

    if norm(y-x,inf) < TOL
        exitflag = 0;
        return
    end

end

exitflag = -1;