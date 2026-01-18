clc;
clear all;

% =========================================================================
% 2D STEADY-STATE HEAT CONDUCTION IN A RECTANGULAR PLATE
% =========================================================================
% This script solves the 2D steady-state heat conduction equation using
% analytical Fourier series solution for a rectangular plate with specified
% boundary conditions.
%
% Boundary Conditions:
%   - Three edges maintained at temperature T1
%   - One edge maintained at temperature T2
%   - No heat generation within the plate
%
% Solution Method:
%   - Separation of variables with Fourier series expansion
%   - Uses hyperbolic sine functions to satisfy boundary conditions
%
% Author: Sofia Vega
% Date: January 2026
% License: MIT
% =========================================================================

%% PLATE GEOMETRY

% Dimensions of the rectangular plate
L = 2;  % Length of the plate [m]
W = 2;  % Width of the plate [m]

%% EVALUATION POINT

% Coordinates where temperature is being evaluated
x_eval = 1;  % x-coordinate of the evaluation point [m]
y_eval = 1;  % y-coordinate of the evaluation point [m]

%% NUMERICAL PARAMETERS

% Number of terms in the Fourier series expansion
% Higher N provides better accuracy but increases computation time
N = 6;  % Number of terms [-]

% Grid resolution for visualization
Nx = 50; % Number of grid points in x-direction
Ny = 50; % Number of grid points in y-direction

%% BOUNDARY CONDITIONS

% Temperature at the plate boundaries
T1 = 280;  % Temperature at three boundaries [K]
T2 = 320;  % Temperature at one boundary (y = W) [K]

%% GRID GENERATION

% Create spatial grid for the plate
x = linspace(0, L, Nx);
y = linspace(0, W, Ny);
[X, Y] = meshgrid(x, y);

%% TEMPERATURE FIELD CALCULATION

% Initialize dimensionless temperature matrix
Theta = zeros(Ny, Nx);

% Compute temperature at each grid point using Fourier series
for i = 1:Ny
    for j = 1:Nx
        theta = 0;
        
        % Sum the Fourier series terms
        for n = 1:N
            % Fourier coefficient (accounts for boundary conditions)
            firstTerm = ((-1)^(n+1) + 1) / n;
            
            % Spatial variation term
            % sin(nπx/L) * sinh(nπy/L) / sinh(nπW/L)
            secondTerm = sin(n * pi * x(j) / L) * ...
                         sinh(n * pi * y(i) / L) / ...
                         sinh(n * pi * W / L);
            
            % Accumulate series sum
            theta = theta + firstTerm * secondTerm;
        end
        
        % Apply normalization factor
        Theta(i, j) = (2 / pi) * theta;
    end
end

% Convert dimensionless temperature to actual temperature
% T(x,y) = θ(x,y) * (T2 - T1) + T1
T = Theta * (T2 - T1) + T1;

%% TEMPERATURE AT SPECIFIC POINT

% Calculate dimensionless temperature at evaluation point
theta_eval = 0;
for n = 1:N
    firstTerm = ((-1)^(n+1) + 1) / n;
    secondTerm = sin(n * pi * x_eval / L) * ...
                 sinh(n * pi * y_eval / L) / ...
                 sinh(n * pi * W / L);
    theta_eval = theta_eval + firstTerm * secondTerm;
end
theta_eval = (2 / pi) * theta_eval;

% Calculate actual temperature at evaluation point
T_eval = theta_eval * (T2 - T1) + T1;

%% DISPLAY RESULTS

fprintf('========================================\n');
fprintf('2D PLATE HEAT CONDUCTION ANALYSIS\n');
fprintf('========================================\n');
fprintf('Plate Dimensions: %.1f m x %.1f m\n', L, W);
fprintf('Boundary Temperatures:\n');
fprintf('  T1 (three edges): %.1f K\n', T1);
fprintf('  T2 (one edge): %.1f K\n', T2);
fprintf('Series Terms: %d\n\n', N);

fprintf('Temperature at point (%.1f, %.1f):\n', x_eval, y_eval);
fprintf('  Dimensionless θ: %.4f\n', theta_eval);
fprintf('  Actual T: %.2f K\n', T_eval);
fprintf('========================================\n\n');

%% VISUALIZATION

figure('Position', [100, 100, 1200, 500]);

% Plot 1: Dimensionless Temperature Distribution
subplot(1, 2, 1);
surf(X, Y, Theta, 'EdgeColor', 'none');
colorbar;
title('Dimensionless Temperature θ(x, y)');
xlabel('x (m)');
ylabel('y (m)');
zlabel('θ');
view(2); % Top view

% Highlight the evaluation point
hold on;
scatter3(x_eval, y_eval, theta_eval, 100, 'r', 'filled');
text(x_eval, y_eval, theta_eval, ...
     ['  θ(', num2str(x_eval), ', ', num2str(y_eval), ') = ', ...
      sprintf('%.2f', theta_eval)], ...
     'Color', 'w', 'FontSize', 12);
hold off;

% Plot 2: Actual Temperature Distribution
subplot(1, 2, 2);
surf(X, Y, T, 'EdgeColor', 'none');
colorbar;
title('Temperature Distribution T(x, y)');
xlabel('x (m)');
ylabel('y (m)');
zlabel('Temperature (K)');
view(2); % Top view

% Highlight the evaluation point
hold on;
scatter3(x_eval, y_eval, T_eval, 100, 'r', 'filled');
text(x_eval, y_eval, T_eval, ...
     ['  T(', num2str(x_eval), ', ', num2str(y_eval), ') = ', ...
      sprintf('%.2f', T_eval), ' K'], ...
     'Color', 'w', 'FontSize', 12);
hold off;
