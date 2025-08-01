% =========================================================================
% FUNCTION NAME: GaussSeidel_IP.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:  This function performs Gauss-Seidel iterations on the
% primal-dual system accounting for sparsity. 
%
% INPUTS:
%   - A: square nxn matrix
%   - b: nx1 rhs vector
%   - x0: initial guess
%   - TOL: tolerance for convergence
%   - MAXIT: maximum number of iterations
%   - p: prediction horizon
%   - m: number of constraints
%
% OUTPUTS:
%   - x: solution to the linear system
%   - k: iterations to converge
%
% EXIT FLAGS:
%   - exitflag  = -1 means algorithm did not converge,
%               = 0 means the algorithm did converge.
% =========================================================================

function [x,k,exitflag] = GaussSeidel_IP(A,b,x0,TOL,MAXIT,p,m)

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
            x(i) = (b(i) - sum1 - sum2)/A(i,i); 
        elseif i <= p+m % Phase 2
            sum1 = 0;
            for j = 1 : p
                sum1 = sum1 + A(i,j)*x(j);
            end
            x(i) = b(i) - sum1; 
        else % Phase 3
            x(i) = (b(i) - A(i,i-m)*x(i-m))/A(i,i); 
        end
    end

    if norm(y-x,inf) < TOL
        exitflag = 0;
        return
    end

end

exitflag = -1;