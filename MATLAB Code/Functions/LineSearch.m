function [alpha] = LineSearch(v, p_v, tau)

alpha = 1e-2; % through trial and error, this is a good value 
while (v + alpha*p_v < (1-tau)*v)
    alpha = 0.95*alpha;
end

