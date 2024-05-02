close all
clear

%% extra parameters
depth = 2700;
h=linspace(0,depth,depth+1);
vAvg = 30; % from literature, m/s
g0 = 9.81; % m/s^2

%% calculate density of water at depth
% density calculations
depthShallow = linspace(0,1000,1000+1);
rhoShallow = 3*10^(-6)*depthShallow+1025; % kg/m^3
rhoDeep = ones(1,depth-1000)*1028; % kg/m^3
rhoWater = [rhoShallow,rhoDeep];

% pressure calculations
pATM = 0.101325; % MPa
pWater = zeros(1,depth+1);
pWater(1) = pATM;
pHydroW = 0;
for i=1:depth
    pHydroWi = rhoWater(i)*g0*10^(-6); % by convension deltaH=1;
    pHydroW = pHydroW + pHydroWi;
    pWater(i+1) = pHydroW;
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
rhoCO2ref = [0 10 30 50 85 110 145 235 500 620 700 790 860 905 935]; %kg/m^3
pCO2ref = [0 1 2 3 4 5 6 7 8 9 10 15 20 25 30]; % MPa

% create fit of p and rho data
[xData, yData] = prepareCurveData( pCO2ref, rhoCO2ref);
ft = 'linearinterp';
opts = fitoptions( 'Method', 'LinearInterpolant' );
opts.ExtrapolationMethod = 'linear';
opts.Normalize = 'on';
[rhoVsPFit,~] = fit( xData, yData, ft, opts );

% test varying values of applied pressure
pApplied = 8; % MPa
pC = pApplied;
pCFxnDepth = zeros(1,depth+1);
pCFxnDepth(1) = pC;
pLoss = 0.5*vAvg*0.016/depth; % MPa/m
pHydroAll=zeros(1,depth);
pLossAll=zeros(1,depth);
rhoAll=zeros(1,depth);

% iteration to calculate CO2 pressure at depth
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

%% graph
figure
plot(pCFxnDepth,-h,'LineWidth',1.5)
xlabel('Pressure of CO2 (MPa)')
ylabel('Depth (m)')
hold on

% constraint: injection pressure
pInjection = 24; % MPa requirement - preliminary finding that needs more research
plot([pInjection,pInjection], [0,-depth],'LineWidth',1.5)

% constraint: maintain supercritical fluid
pSupercritical = 7.37; % MPa requirement
plot([pSupercritical, pSupercritical], [0, -depth],'LineWidth',1.5)

plot(pWater,-h,'LineWidth',1.5)
% contsraint: no pipe explosion
% constraint: no pipe implosion
pDifferential = pWater - pCFxnDepth;
pDiffMax = max(pDifferential);
pDiffMin = min(pDifferential);

legend('CO2 Pressure Chart','Injection MINIMUM Requirment', ...
    'Supercritical MINIMUM Requirement','Water Pressure Chart')

