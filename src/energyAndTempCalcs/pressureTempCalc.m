close all
clear

%% extra parameters
depth = 2700;
h = linspace(0,depth,depth+1);
vAvg = 30; % from literature, m/s
g0 = 9.81; % m/s^2

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
    pWater(i+1) = pWater(i) + rhoWater(i)*g0*10^(-6); % by convention deltaH=1;
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
pLoss = 0.5 * vAvg * 0.016 / depth; % MPa/m
pHydroAll = zeros(1,depth);
pLossAll = zeros(1,depth);
rhoAll = zeros(1,depth);

%% iteration to calculate CO2 pressure at depth
for i=1:depth
    rhoi = rhoVsPFit(pC); % kg/m^3
    rhoAll(i)=rhoi;
    % note: deltahi=1 by convension
    pHydroi = g0*rhoi*10^(-6); % MPa, iterative hydrostatic pressure
    pLossi = pLoss*rhoi*10^(-6); % MPa, iterative pressure loss in pipe
    pC = pC - pLossi + pHydroi;
    pCFxnDepth(i+1) = pC;
    pHydroAll(i)  = pHydroi;
    pLossAll(i) = pLossi;
end

%% Ehina's Code for temperature

%% Temperature Inputs
% CO2_pressure_values= pCFxnDepth;
deltaP = diff(pCFxnDepth);
rhoCO2 = rhoAll;
m = rhoCO2* vAvg;
h_new = 162.3 ; % some constant (for now)
C = 0.709; % some constant (for now)
Toc = 290.15; % some constant (for now)
deltaZ = 1; % try a constant for now
R =  0.009;
k = 45; 
outer_radius_pipe = 0.5; % m
thickness_pipe = 0.01; % m
inner_radius_pipe = outer_radius_pipe - thickness_pipe;
A_pipe = pi*(outer_radius_pipe^2 - inner_radius_pipe); % cross sectional area of pipe


%% Temperature Calculation
rho_pipe = 7900;
c_pipe = 500;
heat_capacity_pipe = A_pipe * deltaZ * rho_pipe * c_pipe;

%create empty lists to populate
TCO2 = zeros(1,depth+1);
TS = zeros(1,depth+1); 

%initialize TCO2 and TS at the surface of the ocean
TCO2(1) = 304.25; %minimum value to maintain supercritical state
TS(1) = 200;

for i = 1:depth
    heat_capacity_CO2 = m(i) * C;
    heat_xfer_from_CO2_to_pipe = 2*pi*R*k*(TCO2(i) - TS(i));

    net_heat_xfer_into_CO2  = -heat_xfer_from_CO2_to_pipe + g0*deltaZ + deltaP(i)/rhoCO2(i);
    net_heat_xfer_into_pipe =  heat_xfer_from_CO2_to_pipe - h_new*A_pipe*(TS(i) - Toc);

    % net heat xfer = heat capacity * ( T(i+1) - T(i) )
    TCO2(i+1) = ( net_heat_xfer_into_CO2  + heat_capacity_CO2  * TCO2(i) ) / heat_capacity_CO2; 
    TS(i+1)   = ( net_heat_xfer_into_pipe + heat_capacity_pipe * TS(i)   ) / heat_capacity_pipe; 
end 


%% Plot temperature and pressure of CO2 as a function of depth
figure
%Pressure
x = h;
y1 = pCFxnDepth;  
y2 = TCO2;

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

plot(h,pWater,'LineWidth',1.5, 'Color','#EDB120', 'LineStyle',"-")
% contsraint: no pipe explosion
% constraint: no pipe implosion
pDifferential = pWater - pCFxnDepth;
pDiffMax = max(pDifferential);
pDiffMin = min(pDifferential);

legend('CO2 Pressure Chart','Injection MINIMUM Requirment', ...
    'Supercritical MINIMUM Requirement','Water Pressure Chart')

