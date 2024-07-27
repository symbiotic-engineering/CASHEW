close all
clear

supercritical = false; % flag for supercritical or liquid CO2

%% parameters
% depth
depth = 2700;
N = 1000; % number of depth steps

% environmental constants
g = 9.81; % m/s^2
P_atm = 101325; % Pa
k_darcy = 10e-6; % seabed permeability, d (darcy) - from Ch3 here: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007GL031560
k = k_darcy * 1e-12; % seabed permeability, m^2

% pipe thermals
h_out    = 50; % convective heat transfer coefficient from pipe to ocean - some constant (for now) [W/(m^2 C)]
h_in     = 50;   % convective heat transfer coefficient from CO2 to pipe [W/(m^2 C)]
T_oc     = 17;    % temperature of ocean - some constant (for now) [C]
k_pipe   = 45;    % thermal conductivity steel [W/(m C)]
k_insu   = 0.015; % thermal conductivity insulation [W/(m C)] - from https://www.aerogel.com/wp-content/uploads/2021/08/Spaceloft-Subsea-Datasheet.pdf
rho_pipe = 7900;  % density steel [kg/m^3]
c_pipe   = 500;   % specific heat steel [J / (kg C)]
if supercritical
    P_heat   = 150;    % heating of CO2 [W] - guess and check until temp graph exceeds supercritical requirement
else
    P_heat = 0;
end

% CO2 properties
P_supercritical = 7.38 * 1e6; % [Pa] requirement
T_supercritical = 31; % minimum temperature of CO2 to maintain supercritical state [C]
if supercritical
    T_CO2_surface = T_supercritical;    
    c_CO2         = 0.709; % specific heat of CO2 - some constant (for now) [J / (kg C)]
else
    % liquid CO2
    T_CO2_surface = 0; % [C] 
    c_CO2 = 2470; % https://www.engineeringtoolbox.com/carbon-dioxide-d_1000.html
end

% pipe dimensions
outer_radius_pipe = 0.14;  % [m]
thickness_pipe    = 0.01; % [m]
thickness_insulation = 0.01; % [m] - from aerogel link above

inner_radius_pipe = outer_radius_pipe - thickness_pipe;
outer_radius_insu = outer_radius_pipe + thickness_insulation;


%% calculate CO2 surface pressure required for different massflows and different injection depths

mdot_Mt_yr = [1 5:5:40]; % massflow of CO2 in megatons/yr
sec_per_yr = 365.25 * 24 * 60 * 60;
mdot = [1 10 25:25:200];%mdot_Mt_yr * 1e9 / sec_per_yr; % convert from megatons/yr to kg/s

injection_depths = 100 : 100 : 400;

P_surface_required = zeros(length(mdot),length(injection_depths));
power = zeros(size(P_surface_required));
for i = 1:length(mdot)
    for j = 1:length(injection_depths)
         [P,~,~,rho]= pressures(depth, P_atm, g, mdot(i), inner_radius_pipe, injection_depths(j), k, N, supercritical);
         P_surface_required(i,j) = P;
         power(i,j) = P / rho(1) * mdot(i);
    end
end


legend_text = cellstr(num2str(injection_depths'));
figure
subplot 121
plot(mdot, power/1e3)
xlabel('Desired Massflow of CO2 kg/s')%(Mt/yr)')
ylabel('Required WEC Power (kW)')
title('Power')
leg = legend(legend_text);
title(leg,'Injection Depth (m)')
improvePlot

subplot 122
%plot(mdot_Mt_yr,P_surface_required/1e6)
plot(mdot,P_surface_required/1e6)
if supercritical
    hold on
    plot([0 max(mdot_Mt_yr)],P_supercritical*[1 1]/1e6,'k--')
    legend_text = [legend_text;'Supercritical'];
end
xlabel('Desired Massflow of CO2 kg/s')%(Mt/yr)')
ylabel('Required Pressure at Surface (MPa)')
leg = legend(legend_text);
title(leg,'Injection Depth (m)')
title('Pressure')
improvePlot

%% single massflow and injection depth to use for rest of analysis
mdot_main = mdot(2);
injection_depth_main = 200;
[P_surface_required,P_bottom,...
 P_vs_depth, rho_vs_depth,P_water,rho_water,...
 h_under_seafloor,deltaZ,z] = pressures(depth, P_atm, g, mdot_main, inner_radius_pipe, injection_depth_main, k, N, supercritical);

%% calculate allowable pressure in a pipe
% contsraint: no pipe explosion
% constraint: no pipe implosion
pDifferential = P_water - P_vs_depth;
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
                    thickness_pipe,thickness_insulation,T_CO2_surface,T_oc,...
                    k_pipe,k_insu,h_in,h_out,P_heat,rho_pipe,c_CO2,c_pipe,...
                    h_under_seafloor,deltaZ,g,N);

%% Plot temperature and pressure and density of CO2 as a function of depth
%Pressure
figure
plot(z,P_vs_depth/1e6,'LineWidth', 1.5,'DisplayName','CO2 Pressure');
hold on
if supercritical
    % constraint: maintain supercritical fluid
    plot([0, depth], [P_supercritical, P_supercritical]/1e6,'--k','LineWidth',1.5,'DisplayName','Supercritical Minimum Requirement')
end
plot(z,P_water/1e6,'LineWidth',1.5,'DisplayName','Water Pressure')

title('Pressure as a Function of Depth')
ylabel('Pressure (MPa)')
xlabel('Depth (m)')
legend
improvePlot

% density
figure
plot(z,rho_vs_depth,'LineWidth', 1.5);
hold on
plot(z,rho_water,'LineWidth',1.5)

title('Density as a Function of Depth')
ylabel('Density (kg/m^3)')
xlabel('Depth (m)')
legend('CO2 Density','Water Density')
improvePlot

% temperature
figure
plot(z,T_CO2, 'LineWidth', 1.5,'DisplayName','CO2 Temperature')
if supercritical
    hold on
    plot([0 max(z)],T_supercritical*[1 1],'DisplayName','Supercritical');
end
ylabel('Temperature (C)')
xlabel('Depth (m)')
legend
improvePlot

%% Plot all temperatures as a function of depth

figure
plot(z, T_CO2, z, T_s_in, '--', z, T_s_out, '-.', z, T_insu,'LineWidth',1.5)
legend_text = {'CO2','Steel pipe inside','Steel pipe outside','Insulation outside'};
if supercritical
    hold on
    plot([0 max(z)],T_supercritical*[1 1],'k--')
    legend_text = [legend_text,'Supercritical'];
end
legend(legend_text)
xlabel('Depth (m)')
ylabel('Temperature (C)')
improvePlot

function [P_surface_required,P_bottom,...
        P_vs_depth, rho_vs_depth,P_water,rho_water,...
        h_under_seafloor,deltaZ,z] = pressures(depth, P_atm, g, mdot, inner_radius_pipe, ...
                                                distance_below_seafloor, k, N, supercritical)

    h_under_seafloor = depth + distance_below_seafloor;
    deltaZ = h_under_seafloor / N;
    z = linspace(0,h_under_seafloor,N);

    % water density change with depth
    water_shallow_threshold = 1000; % m
    depth_shallow = z(z <= water_shallow_threshold); % m
    rho_shallow = 3*10^(-6) * depth_shallow + 1025; % kg/m^3
    depth_deep = z(z > water_shallow_threshold); % m
    rho_deep = ones(size(depth_deep)) * 1028; % kg/m^3
    rho_water = [rho_shallow, rho_deep];
    
    %% calculate water pressure at depth
    P_water = zeros(1,N);
    P_water(1) = P_atm; % Pa
    for i=1:N-1
        P_water(i+1) = P_water(i) + rho_water(i) * g * deltaZ; % by convention deltaH=1;
    end

    D = 2*inner_radius_pipe;
    
    pore_multiplier = 1; % factor of 1 to 1.8 https://petrowiki.spe.org/Methods_to_determine_pore_pressure
    P_pore_bottom = pore_multiplier * P_water(end);
    P_frack = 2 * 1e6 * 0.35 * (distance_below_seafloor/100)^1.5; % rough fit of Fig 3 https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007GL031560
    P_bottom_required = P_pore_bottom + P_frack;

    [P_surface_required,P_bottom,...
        P_vs_depth, rho_vs_depth] = get_P_surface(mdot,depth,h_under_seafloor,P_bottom_required,g,D,k,N,supercritical);
end
