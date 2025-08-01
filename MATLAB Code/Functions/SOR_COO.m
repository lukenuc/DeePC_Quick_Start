
% =========================================================================
% FUNCTION NAME: SOR_COO.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:  Performs SOR where the sparse matrix is represented in COO
% format. 
%
% INPUTS:
%   - arow: row index vector
%   - acol: column index vector
%   - aval: value vector
%   - b: RHS vector
%   - x0: initial guess
%   - omega: relaxation parameter
%   - TOL: tolerance for convergence
%   - MAXIT: maximum number of iterations
%
% OUTPUTS:
%   - x: solution to the linear system
%   - k: iterations to converge
%
% EXIT FLAGS: = -1, did not converge
%             = 0, did converge
% 
% =========================================================================

function [x,k,exitflag] = SOR_COO(arow,acol,aval,b,x0,omega,TOL,MAXIT)

x = x0(:);
n = length(b);
l = length(aval); 

for k = 1 : MAXIT
    y = x; % Save previous iteration
    cnt = 1; 
    for i = 1 : n
        sum1 = 0; sum2 = 0; aii = 0; 
        while (cnt <= l & arow(cnt) == i)
            col = acol(cnt); a = aval(cnt);
            if col < i 
                sum1 = sum1 + a*x(col);
            elseif col > i                 
                sum2 = sum2 + a*y(col); 
            else
                aii = a; 
            end
            cnt = cnt + 1; 
        end
        x(i) = (omega*b(i) - omega*sum1 - omega*sum2 + (1-omega)*aii*y(i))/aii;
    end

    if norm(y-x,inf) < TOL
        exitflag = 0;
        return
    end

end

exitflag = -1;