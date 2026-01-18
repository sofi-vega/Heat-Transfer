clc;
clear;

% =========================================================================
% MULTILAYER CYLINDRICAL HEAT TRANSFER CALCULATOR
% =========================================================================
% This script calculates the heat transfer rate per unit length through a
% multilayer cylindrical pipe with convection at both inner and outer surfaces.
%
% The calculation uses thermal resistance network analysis for:
%   - Inner convection resistance
%   - Conduction resistance through each layer
%   - Outer convection resistance
%
% Author: Sofía Vega
% =========================================================================

%% INPUT PARAMETERS

% Thermal conductivities for each layer [Btu/(h·ft·°F)]
k = [0.43, 0.06];

% Thickness of each layer [ft]
t = [0.007, 0.025];

% Inner diameter of the pipe [ft]
di = 0.1;

% Temperature at inner surface [°F]
T_inner = 121;

% Ambient temperature at outer surface [°F]
T_outer = 18;

% Inner convection heat transfer coefficient [Btu/(h·ft²·°F)]
h_inner = 85.11;

% Outer convection heat transfer coefficient [Btu/(h·ft²·°F)]
h_outer = 12.48;

%% INITIALIZATION

% Number of insulation layers
n = length(k);

% Total thermal resistance per unit length [(h·°F)/Btu]
R_total = 0;

% Current diameter being processed (updates with each layer)
current_diameter = di;

%% THERMAL RESISTANCE CALCULATIONS

% Inner convection resistance per unit length
% R_conv = 1 / (h * π * D)
R_conv_inner = 1 / (h_inner * pi * di);
R_total = R_total + R_conv_inner;

% Conduction resistance through each cylindrical layer
% R_cond = ln(r_outer/r_inner) / (2 * π * k)
for i = 1:n
    % Outer diameter of current layer
    d_outer = current_diameter + 2 * t(i);
    
    % Conduction resistance for this layer per unit length
    R_layer = log(d_outer / current_diameter) / (2 * pi * k(i));
    R_total = R_total + R_layer;
    
    % Update diameter for next layer
    current_diameter = d_outer;
end

% Outer convection resistance per unit length
R_conv_outer = 1 / (h_outer * pi * current_diameter);
R_total = R_total + R_conv_outer;

%% HEAT TRANSFER RATE CALCULATION

% Heat transfer rate per unit length [Btu/(h·ft)]
% Q' = ΔT / R_total
Q_dot = (T_inner - T_outer) / R_total;

%% DISPLAY RESULTS

fprintf('========================================\n');
fprintf('HEAT TRANSFER ANALYSIS RESULTS\n');
fprintf('========================================\n');
fprintf('Heat transfer rate per unit length: %.2f Btu/(h·ft)\n', Q_dot);
fprintf('Total thermal resistance: %.4f (h·°F·ft)/Btu\n', R_total);
fprintf('Final outer diameter: %.4f ft\n', current_diameter);
fprintf('========================================\n');