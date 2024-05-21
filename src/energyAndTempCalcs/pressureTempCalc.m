close all
clear

%% extra parameters
depth = 2700;
g = 9.81; % m/s^2
vel_avg = 3; % m/s

%% calculate density of water at depth
% density calculations
depthShallow = 0:1001; % m
rhoShallow = 3*10^(-6) * depthShallow + 1025; % kg/m^3
rhoDeep = ones(1, depth-1000)*1028; % kg/m^3
rhoWater = [rhoShallow, rhoDeep];

%% water pressure calculations
pATM = 0.101325; % MPa
pWater = zeros(1,depth+1);
pWater(1) = pATM;
for i=1:depth
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

% assume temperature of approx 35 degree C
    % this assumption may be very wrong - we will need to discuss how CO2
    % is maintained at supercritical temperatures during pumping (may
    % because a material requirement for heat conduction/convection
    % coefficient)
rhoCO2ref = [0 10 30 50 85 110 145 235 500 620 700 790 860 905 935]; % kg/m^3
pCO2ref   = [0 1  2  3  4  5   6   7   8   9   10  15  20  25  30 ]; % MPa

% create fit of p and rho data
[xData, yData] = prepareCurveData( pCO2ref, rhoCO2ref);
ft = 'linearinterp';
opts = fitoptions( 'Method', 'LinearInterpolant' );
try
    opts.ExtrapolationMethod = 'linear';
catch
    error('You need R2023 or later to run this code')
end
opts.Normalize = 'on';
[rhoVsPFit,~] = fit( xData, yData, ft, opts );

% test varying values of applied pressure
pApplied = 8; % MPa
pC = pApplied;
pCFxnDepth = zeros(1,depth+1);
pCFxnDepth(1) = pC;
fric_factor = 0.016;
pLoss = 0.5 * vel_avg^2 * fric_factor / depth; % [Pa/(kg/m^3 * m)]
pHydroAll = zeros(1,depth);
pLossAll = zeros(1,depth);
rhoAll = zeros(1,depth);

%% iteration to calculate CO2 pressure at depth
for i=1:depth
    rhoi = rhoVsPFit(pC); % kg/m^3
    rhoAll(i) = rhoi;
    % note: deltahi=1 by convension
    pHydroi = g * rhoi * 10^(-6); % MPa, iterative hydrostatic pressure
    pLossi = pLoss * rhoi * 10^(-6); % MPa, iterative pressure loss in pipe
    pC = pC - pLossi + pHydroi;
    pCFxnDepth(i+1) = pC;
    pHydroAll(i)    = pHydroi;
    pLossAll(i)     = pLossi;
end

%% Ehina's Code for temperature

%% Temperature Inputs
% CO2_pressure_values= pCFxnDepth;
deltaP = diff(pCFxnDepth) * 1e6; % convert MPa to Pa
rhoCO2 = rhoAll;
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
T_CO2   = zeros(1,depth+1); % temperature of CO2

%initialize TCO2 at the surface of the ocean
T_CO2(1) = 31; % minimum value to maintain supercritical state [C]

for i = 1:depth
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
x = 0:depth;
y1 = pCFxnDepth; 
y2 = T_CO2;

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
plot(x,pWater,'LineWidth',1.5, 'Color','#EDB120', 'LineStyle',"-")

% contsraint: no pipe explosion
% constraint: no pipe implosion
pDifferential = pWater - pCFxnDepth;
pDiffMax = max(pDifferential);
pDiffMin = min(pDifferential);

legend('CO2 Pressure Chart','Injection MINIMUM Requirment', ...
    'Supercritical MINIMUM Requirement','Water Pressure')

%% Plot all temperatures as a function of depth

figure
plot(x, T_CO2, x, T_s_in, '--', x, T_s_out, '-.', x, T_insu,'LineWidth',1.5)
legend('CO2','Steel pipe inside','Steel pipe outside','Insulation outside')
xlabel('Depth (m)')
ylabel('Temperature (C)')
