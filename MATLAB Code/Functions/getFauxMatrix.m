
function [H] = getFauxMatrix(N)

m = 6*N; 
H = [rand(N) rand(N,m); rand(m,N) rand(m)];
end

