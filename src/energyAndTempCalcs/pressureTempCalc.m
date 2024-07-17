close all
clear

%% extra parameters
depth = 2700;
g = 9.81; % m/s^2

h_seafloor = depth;
h_under_seafloor = h_seafloor + 200;

%% calculate density of water at depth
% density calculations
depthShallow = 0:1001; % m
rhoShallow = 3*10^(-6) * depthShallow + 1025; % kg/m^3
rhoDeep = ones(1, h_under_seafloor-1000)*1028; % kg/m^3
rhoWater = [rhoShallow, rhoDeep];

%% water pressure calculations
pATM = 0.101325; % MPa
pWater = zeros(1,h_under_seafloor+1);
pWater(1) = pATM;
for i=1:h_under_seafloor-1
    pWater(i+1) = pWater(i) + rhoWater(i) * g * 10^(-6); % by convention deltaH=1;
end

%% calculate allowable pressure in a pipe
% t =;
% r =;
% FOS =; % on stress
%
% % prevent explosion
% sigmaTension =;
% pExpPipe = sigmaTension*t/(r*FOS);
%
% % prevent implosion
% sigmaCompression = ;
% pImpPipe = sigmaCompression*t/(r*FOS);

%% calculate CO2 density at depth
mdot = 1;

[P_surface_required,P_bottom, P_vs_depth, rho_vs_depth] = get_P_surface(mdot,h_seafloor,h_under_seafloor);

% assume temperature of approx 35 degree C
    % this assumption may be very wrong - we will need to discuss how CO2
    % is maintained at supercritical temperatures during pumping (may
    % because a material requirement for heat conduction/convection
    % coefficient)

%% Ehina's Code for temperature



%% Temperature Inputs
% CO2_pressure_values= pCFxnDepth;
deltaP = [0 diff(P_vs_depth)]; % convert MPa to Pa
rhoCO2 = rho_vs_depth;
h_out    = 50; % convective heat transfer coefficient from pipe to ocean - some constant (for now) [W/(m^2 C)]
h_in     = 50;   % convective heat transfer coefficient from CO2 to pipe [W/(m^2 C)]
c_CO2    = 0.709; % specific heat of CO2 - some constant (for now) [J / (kg C)]
T_oc     = 17;    % temperature of ocean - some constant (for now) [C]
deltaZ   = 1;     % distance step [m]
k_pipe   = 45;    % thermal conductivity steel [W/(m C)]
k_insu   = 0.015; % thermal conductivity insulation [W/(m C)] - from https://www.aerogel.com/wp-content/uploads/2021/08/Spaceloft-Subsea-Datasheet.pdf
rho_pipe = 7900;  % density steel [kg/m^3]
c_pipe   = 500;   % specific heat steel [J / (kg C)]
P_heat   = 45;     % heating of CO2 [W]

outer_radius_pipe = 0.5;  % [m]
thickness_pipe    = 0.01; % [m]
thickness_insulation = 0.01; % [m] - from aerogel link above

inner_radius_pipe = outer_radius_pipe - thickness_pipe;
outer_radius_insu = outer_radius_pipe + thickness_insulation;

dA_lat_pipe_outer = 2 * pi * outer_radius_pipe * deltaZ; % lateral area of a sliver of outside of pipe with length deltaZ, m^2
dA_lat_pipe_inner = 2 * pi * inner_radius_pipe * deltaZ; % lateral area of a sliver of inside  of pipe with length deltaZ, m^2
dA_lat_insu_outer = 2 * pi * outer_radius_insu * deltaZ; % lateral area of a sliver of outside of insulation with length deltaZ, m^2

A_cross_pipe = pi * (outer_radius_pipe^2 - inner_radius_pipe); % cross sectional area of pipe, m^2
A_cross_CO2 = pi * inner_radius_pipe^2;                        % cross sectional area of CO2,  m^2

heat_capacity_CO2  = A_cross_CO2  * deltaZ * rhoCO2   * c_CO2;  % mc for a sliver of CO2  with length deltaZ, J/C
heat_capacity_pipe = A_cross_pipe * deltaZ * rho_pipe * c_pipe; % mc for a sliver of pipe with length deltaZ, J/C

conductance_pipe       = dA_lat_pipe_inner * k_pipe / thickness_pipe;       % thermal conductance from pipe inner to pipe outer wall
conductance_insulation = dA_lat_insu_outer * k_insu / thickness_insulation; % thermal conductance from insulation inner to outer wall

therm_resistance_cond_pipe = 1/conductance_pipe;
therm_resistance_cond_insu = 1/conductance_insulation;
therm_resistance_conv_in  = 1/(h_in  * dA_lat_pipe_inner);
therm_resistance_conv_out = 1/(h_out * dA_lat_pipe_outer);
total_thermal_resistance = therm_resistance_cond_pipe + therm_resistance_cond_insu + therm_resistance_conv_in + therm_resistance_conv_out;

%% Temperature Calculation
% create empty lists to populate
T_CO2   = zeros(1,h_under_seafloor+1); % temperature of CO2

%initialize TCO2 at the surface of the ocean
T_CO2(1) = 31; % minimum value to maintain supercritical state [C]

for i = 1:h_under_seafloor
    heat_xfer_from_CO2_to_ocean = 1/total_thermal_resistance * (T_CO2(i) - T_oc);
    net_heat_xfer_into_CO2  = -heat_xfer_from_CO2_to_ocean + g * deltaZ + deltaP(i) / rhoCO2(i) + P_heat;

    % net heat xfer = heat capacity * ( T(i+1) - T(i) )
    T_CO2(i+1) = ( net_heat_xfer_into_CO2  + heat_capacity_CO2(i) * T_CO2(i) ) / heat_capacity_CO2(i);
end

% thermal resistance ratios to find temperature at different points radially
R_ratio_T_s_in  = (therm_resistance_cond_pipe + therm_resistance_cond_insu + therm_resistance_conv_out)/total_thermal_resistance;
R_ratio_T_s_out = (therm_resistance_cond_insu + therm_resistance_conv_out)/total_thermal_resistance;
R_ratio_T_insu  = (therm_resistance_conv_out)/total_thermal_resistance;

T_s_in  = T_oc + (T_CO2 - T_oc) * R_ratio_T_s_in;
T_s_out = T_oc + (T_CO2 - T_oc) * R_ratio_T_s_out;
T_insu  = T_oc + (T_CO2 - T_oc) * R_ratio_T_insu;

%% Plot temperature and pressure of CO2 as a function of depth
figure
%Pressure
x = 1:h_under_seafloor;
y1 = P_vs_depth*1e6; 
y2 = T_CO2(2:end);

yyaxis left
plot(x,y1,'LineWidth', 1.5, 'Color','#0072BD','LineStyle',"-");
title('Pressure and Temperature as a Function of Depth')
ylabel('Pressure of CO2 (MPa)')
xlabel('Depth (m)')

yyaxis right
plot(x,y2, 'LineWidth', 1.5, 'Color','#8040E6','LineStyle',"-")
ylabel('Temperature of CO2 (K)')

yyaxis left
hold on
% constraint: injection pressure
pInjection = 24; % MPa requirement - preliminary finding that needs more research
plot([0,depth],[pInjection,pInjection] ,'LineWidth',1.5, 'Color', '#77AC30', 'LineStyle',"-")

% constraint: maintain supercritical fluid
pSupercritical = 7.37; % MPa requirement
plot( [0, depth], [pSupercritical, pSupercritical],'LineWidth',1.5, 'Color', '#77AC30','LineStyle',"-")
plot(x,pWater(2:end),'LineWidth',1.5, 'Color','#EDB120', 'LineStyle',"-")

% contsraint: no pipe explosion
% constraint: no pipe implosion
pDifferential = pWater(2:end) - P_vs_depth*1e6;
pDiffMax = max(pDifferential);
pDiffMin = min(pDifferential);

legend('CO2 Pressure Chart','Injection MINIMUM Requirment', ...
    'Supercritical MINIMUM Requirement','Water Pressure')

%% Plot all temperatures as a function of depth

figure
plot(x, T_CO2(2:end), x, T_s_in(2:end), '--', x, T_s_out(2:end), '-.', x, T_insu(2:end),'LineWidth',1.5)
legend('CO2','Steel pipe inside','Steel pipe outside','Insulation outside')
xlabel('Depth (m)')
ylabel('Temperature (C)')

function [P_surface_required,P_bottom, P_vs_depth, rho_vs_depth] = get_P_surface(mdot,h_seafloor,h_under_seafloor)
    
    P_surface_guess = 8e6; % 8 MPa
    g = 9.8;
    rho_w = 1000;

    % iterate until P_surface solves the nonlinear equation
    P_bottom_required = rho_w * g * h_under_seafloor;
    P_bottom_error_fcn = @(P_surface_guess) (pressure_vs_depth_fcn(mdot,P_surface_guess,h_seafloor,h_under_seafloor) - P_bottom_required);
    P_surface_required = fzero(P_bottom_error_fcn, P_surface_guess);

    % rerun result with converged P_surface
    [P_bottom, P_vs_depth, rho_vs_depth] = pressure_vs_depth_fcn(mdot,P_surface_required,h_seafloor,h_under_seafloor);

end

function [P_bottom, P_vs_depth, rho_vs_depth] = pressure_vs_depth_fcn(mdot,P_surface,h_seafloor,h_under_seafloor)
    rhoCO2ref = [0 10 30 50 85 110 145 235 500 620 700 790 860 905 935]; % kg/m^3
    pCO2ref   = [0 1  2  3  4  5   6   7   8   9   10  15  20  25  30 ]*1e6; % Pa

    g = 9.8;
    D = 1;

    P_vs_depth = zeros(1,h_under_seafloor);
    P_vs_depth(1) = P_surface;
    rho_vs_depth = zeros(1,h_under_seafloor);
    
    %% iteration to calculate CO2 pressure at depth
    for i=2:h_under_seafloor
        rho_i = interp1(pCO2ref, rhoCO2ref, P_vs_depth(i-1),'linear','extrap'); % density of CO2, kg/m^3
        mu_i = 1; % dynamic viscosity of CO2

        % note: delta_h_i=1 by convention
        deltaP_hydro_i = g * rho_i; % Pa, iterative hydrostatic pressure
        if i < h_seafloor % pipe loss
            deltaP_loss_i = 128 * mdot / (pi * D^4) * mu_i / rho_i^2; % assumes laminar flow so f = 64/Re
        else % injection loss
            fracking = false;
            if fracking
                delta_P_loss_i = P_fracking;
            else % linear Darcy loss based on rock permeability
                C = 1;
                delta_P_loss_i = C * mdot;
            end
        end
        deltaP_i = deltaP_hydro_i - deltaP_loss_i;

        rho_vs_depth(i) = rho_i;
        P_vs_depth(i) = P_vs_depth(i-1) + deltaP_i;
    end
    P_bottom = P_vs_depth(end);

    rho_vs_depth(1) = rho_vs_depth(2);

end

