clc;
clear all;

% =========================================================================
% FIN HEAT TRANSFER ANALYSIS
% =========================================================================
% This script analyzes heat transfer through cylindrical fins made of
% different materials (copper, aluminum, and steel). It calculates and
% plots temperature profiles along the fin length and determines the heat
% transfer rate for each material.
%
% The analysis assumes:
%   - One-dimensional heat conduction along the fin
%   - Constant cross-sectional area
%   - Uniform convection coefficient
%   - Negligible heat loss from the fin tip
%   - Steady-state conditions
%
% Author: Sofia Vega
% Date: January 2026
% License: MIT
% =========================================================================

%% INPUT PARAMETERS

% Fin geometry
D = 0.005; % Diameter of the fin [m]
A = pi * D^2 / 4; % Cross-sectional area of the fin [m²]

% Thermal boundary conditions
Tb = 100; % Base temperature [°C]
Tm = 25; % Ambient temperature [°C]
h = 100; % Convective heat transfer coefficient [W/(m²·K)]

% Thermal conductivity of different materials [W/(m·K)]
k_c = 398; % Copper
k_al = 237; % Aluminum
k_s = 50; % Steel

% Array of thermal conductivities for analysis
k = [k_c, k_al, k_s];

%% PARAMETER CALCULATIONS

% Calculate the fin parameter m for each material
% m = sqrt(hP / kA) where P is the perimeter
% For a cylindrical fin: m = sqrt(h / (k * π * D * A))
m = sqrt(h ./ (k * (pi * D * A)));

% Define the spatial domain along the fin length
x = linspace(0, 0.01, 100); % Length of the fin [m]

%% TEMPERATURE PROFILE CALCULATION

% Preallocate matrix for temperature profiles
Theta = zeros(length(x), length(k));

% Calculate temperature distribution for each material
% Using the equation: θ(x) = (Tb - Tm) * exp(-m*x)
% where θ(x) is the excess temperature at position x
for i = 1:length(k)
    Theta(:, i) = (Tb - Tm) * exp(-m(i) * x);
end

%% VISUALIZATION

% Create figure for temperature profiles
figure;
hold on;

% Define plot aesthetics
colors = ['r', 'g', 'b']; % Red, Green, Blue
labels = {'Copper', 'Aluminum', 'Steel'}; % Material names

% Plot temperature profile for each material
for i = 1:length(k)
    % Add ambient temperature to get actual temperature values
    plot(x, Theta(:, i) + Tm, colors(i), 'LineWidth', 2);
    % Prepare legend information with material properties
    legendInfo{i} = ['Material: ' labels{i} ' (k = ' num2str(k(i)) ' W/(m·K))'];
end

% Enhance graph appearance
title('Temperature Profile along the Fin Length');
xlabel('Fin length (m)');
ylabel('Temperature (°C)');
legend(legendInfo, 'Location', 'best');
grid on;
hold off;

%% HEAT TRANSFER RATE CALCULATION

% Calculate heat transfer rate for each material [W]
% Using the equation: Q = sqrt(h * P * k * A) * (Tb - Tm)
% where P = π * D is the perimeter
Q = sqrt(h * pi * D * k * A) * (Tb - Tm);

%% DISPLAY RESULTS

fprintf('========================================\n');
fprintf('FIN HEAT TRANSFER ANALYSIS RESULTS\n');
fprintf('========================================\n');
fprintf('Fin Diameter: %.4f m\n', D);
fprintf('Base Temperature: %.1f °C\n', Tb);
fprintf('Ambient Temperature: %.1f °C\n', Tm);
fprintf('Convection Coefficient: %.1f W/(m²·K)\n\n', h);

fprintf('Heat Transfer Rate for each material:\n');
for i = 1:length(k)
    fprintf('  %s (k = %.0f W/(m·K)): %.4f W\n', labels{i}, k(i), Q(i));
end
fprintf('========================================\n');
