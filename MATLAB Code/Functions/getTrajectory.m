
% This function takes an input trajectory, and given an LTI system, returns
% the corresponding output trajectory

function [y_traj] = getTrajectory(fcn_dynamics, fcn_meas, u_traj, x0)
    
len = length(u_traj); 
y_traj = zeros(len, 1); 
x = x0; 

for i = 1:len
    u = u_traj(i); 
    y_traj(i,:) = fcn_meas(x, u); 
    x = fcn_dynamics(x, u); 
end

end

