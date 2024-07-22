close all
clear

%% parameters
% depth
depth = 2700;
distance_below_seafloor = 200;
h_under_seafloor = depth + distance_below_seafloor;

% environmental constants
g = 9.81; % m/s^2
P_atm = 101325; % Pa
k_darcy = 10e-6; % seabed permeability, d (darcy) - from Ch3 here: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007GL031560
k = k_darcy * 1e-12; % seabed permeability, m^2

% pipe thermals
h_out    = 50; % convective heat transfer coefficient from pipe to ocean - some constant (for now) [W/(m^2 C)]
h_in     = 50;   % convective heat transfer coefficient from CO2 to pipe [W/(m^2 C)]
T_oc     = 17;    % temperature of ocean - some constant (for now) [C]
deltaZ   = 1;     % distance step [m]
k_pipe   = 45;    % thermal conductivity steel [W/(m C)]
k_insu   = 0.015; % thermal conductivity insulation [W/(m C)] - from https://www.aerogel.com/wp-content/uploads/2021/08/Spaceloft-Subsea-Datasheet.pdf
rho_pipe = 7900;  % density steel [kg/m^3]
c_pipe   = 500;   % specific heat steel [J / (kg C)]
P_heat   = 45;    % heating of CO2 [W]

% supercritical CO2
T_supercritical = 31;    % minimum temperature of CO2 to maintain supercritical state [C]
c_CO2           = 0.709; % specific heat of CO2 - some constant (for now) [J / (kg C)]
P_supercritical = 7.37 * 1e6; % [Pa] requirement

% pipe dimensions
outer_radius_pipe = 0.5;  % [m]
thickness_pipe    = 0.01; % [m]
thickness_insulation = 0.01; % [m] - from aerogel link above

inner_radius_pipe = outer_radius_pipe - thickness_pipe;
outer_radius_insu = outer_radius_pipe + thickness_insulation;


%% calculate water density at depth
% density calculations
depthShallow = 0:1001; % m
rhoShallow = 3*10^(-6) * depthShallow + 1025; % kg/m^3
rhoDeep = ones(1, h_under_seafloor-1000)*1028; % kg/m^3
rhoWater = [rhoShallow, rhoDeep];

%% calculate water pressure at depth
P_water = zeros(1,h_under_seafloor+1);
P_water(1) = P_atm; % Pa
for i=1:h_under_seafloor
    P_water(i+1) = P_water(i) + rhoWater(i) * g; % by convention deltaH=1;
end

%% calculate CO2 density at depth for different massflows

mdot_Mt_yr = [1 20:20:180]; % massflow of CO2 in megatons/yr
sec_per_yr = 365.25 * 24 * 60 * 60;
mdot = mdot_Mt_yr * 1e9 / sec_per_yr; % convert from megatons/yr to kg/s
D = 2*inner_radius_pipe;

P_bottom_required = P_water(end);
P_surface_required = zeros(size(mdot));
for i=1:length(mdot)
    [P_surface_required(i),P_bottom,...
        P_vs_depth, rho_vs_depth] = get_P_surface(mdot(i),depth,h_under_seafloor,P_bottom_required,g,D,k);
end

figure
plot(mdot_Mt_yr,P_surface_required/1e6)
hold on
plot([0 max(mdot_Mt_yr)],P_supercritical*[1 1]/1e6)
xlabel('Desired Massflow of CO2 (Mt/yr)')
ylabel('Required Pressure at Surface (MPa)')
legend('Calculation from density and loss','Supercritical requirement')
improvePlot

%% calculate allowable pressure in a pipe
% contsraint: no pipe explosion
% constraint: no pipe implosion
pDifferential = P_water(2:end) - P_vs_depth;
pDiffMax = max(pDifferential);
pDiffMin = min(pDifferential);
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

%% calculate CO2 and pipe temperature at depth
[T_CO2, T_s_in, T_s_out, T_insu] = temp_model(P_vs_depth,rho_vs_depth,...
                    outer_radius_pipe,inner_radius_pipe,outer_radius_insu,...
                    thickness_pipe,thickness_insulation,T_supercritical,T_oc,...
                    k_pipe,k_insu,h_in,h_out,P_heat,rho_pipe,c_CO2,c_pipe,...
                    h_under_seafloor,deltaZ,g);

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
plot( [0, depth], [P_supercritical, P_supercritical]/1e6,'LineWidth',1.5, 'Color', '#77AC30','LineStyle',"-")
plot(x,P_water(2:end)/1e6,'LineWidth',1.5, 'Color','#EDB120', 'LineStyle',"-")

legend('CO2 Pressure Chart','Injection MINIMUM Requirment', ...
    'Supercritical MINIMUM Requirement','Water Pressure')

%% Plot all temperatures as a function of depth

figure
plot(x, T_CO2(2:end), x, T_s_in(2:end), '--', x, T_s_out(2:end), '-.', x, T_insu(2:end),'LineWidth',1.5)
legend('CO2','Steel pipe inside','Steel pipe outside','Insulation outside')
xlabel('Depth (m)')
ylabel('Temperature (C)')

