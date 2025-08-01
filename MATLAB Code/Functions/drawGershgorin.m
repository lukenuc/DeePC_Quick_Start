
% =========================================================================
% FUNCTION NAME: drawGershgorin.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:   Plots the Gershgorin circles for an input matrix, with
% real values on the horizontal axis and imaginary values on the vertical
% axis. The scale of these axes adjusts to the proper size to fit all
% Gershgorin circles into view. 
%
% INPUTS:
%   - A: square nxn matrix
%
% OUTPUTS: (none)
%
% EXIT FLAGS: exitflag: = -1 means A is not square
%                       = 0 means the function was successful
% 
% =========================================================================

function [exitflag] = drawGershgorin(A)

if (size(A,1) ~= size(A,2))
    exitflag = -1; 
    return
end

n = size(A,1); 
centers = zeros(2, n); 
radii = zeros(1, n); 
colors = rand(n, 3); % random colors
for i = 1:n
    centers(1,i) = real(A(i,i)); % real component
    centers(2,i) = imag(A(i,i)); % imaginary component
    for j = 1:n
        if j ~= i
            radii(i) = radii(i) + abs(A(i,j)); 
        end
    end
end

f = figure;
f.Position = [100 100 500 500]; 
theta = linspace(0, 2*pi, 100); 
x = repmat(centers(1,:)', 1, length(theta)) + radii'.*(cos(theta)); 
y = repmat(centers(2,:)', 1, length(theta)) + radii'.*(sin(theta)); 
x_max = max(max(abs(x))); 
y_max = max(max(abs(x))); 
bound = 1.25*max(x_max, y_max); 

hold on
plot([-bound bound], [0 0], 'k', 'LineWidth', 2)
plot([0 0], [-bound bound], 'k', 'LineWidth', 2)
xlim([-bound bound])
ylim([-bound bound])
for i = 1:n
    fill(x(i,:), y(i,:), colors(i,:), 'FaceAlpha', 0.05);
end
xlabel('Real', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Imaginary', 'FontSize', 14, 'FontWeight', 'bold')
grid minor
hold off

exitflag = 0; 
end

