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
N = 40;                 % [-]   Number of pipe segments
