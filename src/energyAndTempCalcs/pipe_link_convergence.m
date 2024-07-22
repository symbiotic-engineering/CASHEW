clear;clc;close all;
% Initial Conditions
P = 8;                  % [MPa] Intake CO2 Pressure
T = 305;                % [K]   Temperature of CO2

% Environmental Parameters
H_ocean = 2700;         % [m]   Depth of Ocean
H_ground = 200;         % [m]   Distance Underground
g = 9.8;                % [m/s^2]   Gravitational Acceleration
rho_w = 1025;           % [kg/m^3]  Density of Seawater
P_floor = rho_w*g*H_ocean;  % [Pa]  Pressure to overcome at ocean floor
P_floor = P_floor*1e-6  % [MPa] Pressure to overcome at ocean floor

% Piston Parameters
v_amp = 1;              % [m/s] Velocity Amplitude
piston_area = 0.26;     % [m^2] Area of Piston
omega = 1;              % [rad/s]   Wave Frequency

% Pipe Parameters
H = H_ocean + H_ground; % [m]   Total length of pipe
d_pipe = 0.26;          % [m]   Dimater of pipe
A_pipe = d_pipe^2*pi/4; % [m^2] Area of pipe
N = 1:200;            % [-]   Number of pipe segments

mdot = 31710/2700;
P_top = zeros(size(N));

for ii = N
    disp(ii)
    P_top(ii) = get_P_top(P_floor,H,N(ii),g,d_pipe,mdot);
end

plot(N,P_top,'linewidth',2)
xlabel('Number of links')
ylabel('Pressure at top [MPa]')
function P = get_P_top(P,H,N,g,D,mdot)
    for i = 1:N
        P = get_P_above(P,H/N,g,D,mdot);
    end
end

function P_above = get_P_above(P,dH,g,D,mdot)
    rhoCO2ref = [0 10 30 50 85 110 145 235 500 620 700 790 860 905 935]; % kg/m^3
    pCO2ref   = [0 1  2  3  4  5   6   7   8   9   10  15  20  25  30 ]; % MPa
    rho = interp1(pCO2ref, rhoCO2ref, P,'linear','extrap'); % density of CO2, kg/m^3

    Q = mdot / rho; % volume flow rate
    A = pi/4 * D^2;   % area
    v = Q/A;          % velocity

    dP_hydro = g * rho * dH*1e-6; % change in pressure due to gravity
    f = 0.015; % turbulent friction factor, taken for roughness ~5e-4 at high Re on Moody chart.
    dP_fric = 1/2 * rho * v^2 * f/D*dH*1e-6;

    P_above = P - dP_hydro + dP_fric;
end