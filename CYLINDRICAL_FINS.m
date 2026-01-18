clc;
clear;

% =========================================================================
% CYLINDRICAL FINS HEAT TRANSFER ENHANCEMENT ANALYSIS
% =========================================================================
% This script calculates the heat transfer enhancement achieved by adding
% radial fins to a cylindrical surface. It compares the heat dissipation
% with and without fins to quantify the improvement.
%
% The analysis includes:
%   - Fin efficiency calculations for annular fins
%   - Corrected fin dimensions accounting for thickness
%   - Heat transfer comparison between finned and unfinned surfaces
%
% Author: Sofia Vega
% Date: January 2026
% License: MIT
% =========================================================================

%% MATERIAL AND THERMAL PROPERTIES

% Thermal conductivity of fin material [W/(m·K)]
k = 186;

% Convective heat transfer coefficient [W/(m²·K)]
h = 50;

%% GEOMETRIC PARAMETERS

% Cylinder dimensions
L = 0.15;       % Height of the cylinder [m]
D = 0.05;       % Diameter of the cylinder [m]
r = D/2;        % Radius of the cylinder [m]

% Fin dimensions
t = 0.006;      % Thickness of the fins [m]
l = 0.02;       % Radial length of the fins [m]
n = 5;          % Number of fins [-]

%% BOUNDARY CONDITIONS

tm = 300;       % Ambient temperature [K]
tb = 500;       % Base temperature of the cylinder [K]
theta_b = tb - tm; % Temperature difference [K]

%% CORRECTED FIN DIMENSIONS

% Account for fin thickness in calculations
r_ext = r + l;              % Outer radius of fins [m]
r_ext_c = r_ext + t/2;      % Corrected outer radius including thickness [m]
lc = l + t/2;               % Corrected radial length [m]

%% AREA CALCULATIONS

% Fin profile area (cross-sectional area for heat conduction)
Ap = lc * t;                % [m²]

% Total heat transfer surface area of all fins
% Includes both top and bottom surfaces of annular fins
Af = n * pi * (r_ext_c^2 - r^2) * 2; % [m²]

% Cylinder base area without fins
A_base = 2 * pi * r * L;    % [m²]

% Cylinder area adjusted for fin locations
A_cylinder_without_fins = 2 * pi * r * (L - n * t); % [m²]

%% FIN EFFICIENCY CALCULATION

% Radius ratio for annular fin efficiency correlation
ratio = r_ext_c / r;

% Fin parameter for efficiency calculation
% ξ = sqrt(h / (k * Ap)) * lc^(3/2)
cita = sqrt(h / (k * Ap)) * lc^(3/2);

% Fin efficiency (empirical value or from chart)
% For this geometry, assumed efficiency
eta = 0.95;                 % [-]

%% HEAT TRANSFER CALCULATIONS

% Heat transfer WITHOUT fins [W]
% Q = h * A * ΔT
Q_nofin = h * A_base * theta_b;

% Heat transfer WITH fins [W]
% Q = h * (A_unfinned + η * A_fins) * ΔT
Q_fin = h * (A_cylinder_without_fins + eta * Af) * theta_b;

% Heat transfer enhancement due to fins [W]
delta_Q = Q_fin - Q_nofin;

% Heat transfer improvement percentage
improvement_percent = (delta_Q / Q_nofin) * 100;

%% DISPLAY RESULTS

fprintf('========================================\n');
fprintf('CYLINDRICAL FIN ANALYSIS RESULTS\n');
fprintf('========================================\n');
fprintf('Cylinder Dimensions:\n');
fprintf('  Height: %.3f m\n', L);
fprintf('  Diameter: %.3f m\n', D);
fprintf('  Number of fins: %d\n\n', n);

fprintf('Fin Geometry:\n');
fprintf('  Thickness: %.4f m\n', t);
fprintf('  Radial length: %.4f m\n', l);
fprintf('  Fin efficiency: %.2f%%\n\n', eta * 100);

fprintf('Thermal Performance:\n');
fprintf('  Heat transfer WITHOUT fins: %.2f W\n', Q_nofin);
fprintf('  Heat transfer WITH fins: %.2f W\n', Q_fin);
fprintf('  Increment in heat transfer: %.2f W\n', delta_Q);
fprintf('  Improvement: %.1f%%\n', improvement_percent);
fprintf('========================================\n');
