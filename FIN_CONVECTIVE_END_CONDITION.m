clc;
clear all;

% =========================================================================
% FIN HEAT TRANSFER ANALYSIS - CONVECTIVE END CONDITION
% =========================================================================
% This script analyzes heat transfer through a fin with uniform cross-section
% considering convection at the fin tip (finite length with convective end).
%
% Boundary Conditions:
%   - Base: Fixed temperature Tb
%   - Tip: Convective heat transfer to ambient
%   - Lateral surface: Convection along entire length
%
% Solution includes:
%   - Exact analytical solution for temperature distribution
%   - Harper-Brown corrected length approximation
%   - Heat flux calculations for both methods
%
% Author: Sofia Vega
% Date: January 2026
% License: MIT
% =========================================================================

%% THERMAL PARAMETERS

% Boundary temperatures
Tb = 100;       % Temperature at the base of the fin [°C]
T_inf = 20;     % Ambient temperature [°C]

% Material and heat transfer properties
k = 237;        % Thermal conductivity of the material [W/(m·K)]
h = 35;         % Convection coefficient [W/(m²·K)]

% Fin length
L = 0.03;       % Length of the fin [m]

%% FIN GEOMETRY DEFINITION

% OPTION 1: Rectangular Fin (commented out)
% w = 7;                      % Width of the fin [m]
% t = 0.1;                    % Thickness of the fin [m]
% P = 2*w + 2*t;              % Perimeter [m]
% A_cond = w*t;               % Cross-sectional area [m²]
% m_aprox = sqrt((2*h)/(k*t));% Approximate m for thin rectangular fin
% Lc = L + (t/2);             % Harper-Brown corrected length [m]

% OPTION 2: Cylindrical Fin (active)
r = 1.25e-3;                % Radius of the fin [m]
P = 2*pi*r;                 % Perimeter of the cylindrical fin [m]
A_cond = pi*r^2;            % Cross-sectional area [m²]
Lc = L + (r/2);             % Harper-Brown corrected length [m]

%% FIN PARAMETER CALCULATION

% Calculate fin parameter m
% m = sqrt(hP / kA)
m = sqrt((h * P) / (A_cond * k));
m_values = [m];             % Array for plotting multiple m values

%% SPATIAL DISCRETIZATION

% Create spatial grid for temperature profile
x = linspace(0, L, 100);    % Position along fin length [m]

% Initialize temperature storage matrix
T_values = zeros(length(m_values), length(x));

%% HEAT FLUX CALCULATIONS

% Exact solution for heat flux with convective end condition
% Q = kA * (Tb - T∞) * m * [tanh(mL) + h/(km)] / [1 + h*tanh(mL)/(km)]
numerator = tanh(m*L) + (h/(k*m));
denominator = 1 + (h*tanh(m*L)/(k*m));
Q = k * A_cond * (Tb - T_inf) * m * (numerator / denominator);

% Harper-Brown approximation (insulated tip with corrected length)
% Q ≈ kA * (Tb - T∞) * m * tanh(m*Lc)
Q_aprox = k * A_cond * (Tb - T_inf) * m * tanh(m*Lc);

%% TEMPERATURE PROFILE CALCULATION

% Calculate temperature distribution along the fin
for i = 1:length(m_values)
    m_current = m_values(i);
    
    % Analytical solution for temperature profile
    % T(x) = (Tb - T∞) * [cosh(m(L-x)) + (h/(km))sinh(m(L-x))] / 
    %                     [cosh(mL) + (h/(km))sinh(mL)] + T∞
    numerator_T = cosh(m_current*(L - x)) + ...
                  (h/(k*m_current)) * sinh(m_current*(L - x));
    denominator_T = cosh(m_current*L) + ...
                    (h/(k*m_current)) * sinh(m_current*L);
    
    T = (Tb - T_inf) * (numerator_T / denominator_T) + T_inf;
    T_values(i, :) = T;
end

%% DISPLAY RESULTS

fprintf('========================================\n');
fprintf('FIN ANALYSIS - CONVECTIVE END CONDITION\n');
fprintf('========================================\n');
fprintf('Fin Geometry:\n');
fprintf('  Type: Cylindrical\n');
fprintf('  Radius: %.4f mm\n', r*1000);
fprintf('  Length: %.3f m\n', L);
fprintf('  Corrected length (Lc): %.4f m\n\n', Lc);

fprintf('Thermal Properties:\n');
fprintf('  Base temperature: %.1f °C\n', Tb);
fprintf('  Ambient temperature: %.1f °C\n', T_inf);
fprintf('  Thermal conductivity: %.1f W/(m·K)\n', k);
fprintf('  Convection coefficient: %.1f W/(m²·K)\n\n', h);

fprintf('Fin Parameter:\n');
fprintf('  m = %.4f m⁻¹\n\n', m);

fprintf('Heat Transfer Results:\n');
fprintf('  Exact solution: %.6f W\n', Q);
fprintf('  Harper-Brown approximation: %.6f W\n', Q_aprox);
fprintf('  Relative error: %.2f%%\n', abs((Q - Q_aprox)/Q)*100);
fprintf('========================================\n\n');

%% VISUALIZATION

figure('Position', [100, 100, 800, 600]);
hold on;

% Plot temperature profile for each m value
for i = 1:length(m_values)
    plot(x, T_values(i, :), 'LineWidth', 2, ...
         'DisplayName', ['m = ' num2str(m_values(i), '%.2f') ' m⁻¹']);
end

% Format plot
xlabel('Position along fin, x (m)', 'FontSize', 12);
ylabel('Temperature, T(x) (°C)', 'FontSize', 12);
title('Temperature Profile along Fin - Case 3: Convective End', 'FontSize', 14);
legend('show', 'Location', 'best');
grid on;
hold off;

% Add annotations for boundary conditions
ylim_curr = ylim;
text(0, Tb, sprintf('  Base: %.0f°C', Tb), ...
     'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'r');
text(L, T_values(1,end), sprintf('  Tip (convective)'), ...
     'VerticalAlignment', 'top', 'FontSize', 10, 'Color', 'b');
