clear;clc;close all;
% Environmental Parameters
H_ocean = 2700;         % [m]   Depth of Ocean
H_ground = 200;         % [m]   Distance Underground
g = 9.8;                % [m/s^2]   Gravitational Acceleration
rho_w = 1025;           % [kg/m^3]  Density of Seawater
P_ocean = rho_w*g*H_ocean;  % [Pa]  Pressure at the ocean floor
P_frac = 2;             % [MPa] Pressure required to break rocks
P_floor = P_ocean*1e-6 + P_frac;    % [MPa] Total pressure to frac

% Pipe Parameters
H = H_ocean + H_ground; % [m]   Total length of pipe
d_pipe = 0.26;          % [m]   Dimater of pipe
A_pipe = d_pipe^2*pi/4; % [m^2] Area of pipe
N = [40,2700];            % [-]   Number of pipe segments

mdot = 300;
P_top = zeros(size(N));

for ii = 1:length(N)
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
    pCO2ref = [3.952702582,5.000000552,6.013513969,6.993242831,8.006756248,9.020269664,9.999999816,19.99999834,29.99999945];
    rhoCO2ref = [933.1932527,943.6974547,950.0000401,958.4033536,962.6050504,968.9075556,977.3109493,1021.428566,1055.04202];
    rho = interp1(pCO2ref, rhoCO2ref, P,'linear','extrap'); % density of CO2, kg/m^3

    Q = mdot / rho; % volume flow rate
    A = pi/4 * D^2;   % area
    v = Q/A;          % velocity

    dP_hydro = g * rho * dH*1e-6; % change in pressure due to gravity
    f = 0.015; % turbulent friction factor, taken for roughness ~5e-4 at high Re on Moody chart.
    dP_fric = 1/2 * rho * v^2 * f/D*dH*1e-6;

    P_above = P - dP_hydro + dP_fric;
end