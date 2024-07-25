clear;clc;close all;

% Initial Conditions
P = 8;                  % [MPa] Intake CO2 Pressure
T = 273.15;                % [K]   Temperature of CO2

% Environmental Parameters
H_ocean = 2700;         % [m]   Depth of Ocean
H_ground = 200;         % [m]   Distance Underground
g = 9.8;                % [m/s^2]   Gravitational Acceleration
rho_w = 1025;           % [kg/m^3]  Density of Seawater
P_ocean = rho_w*g*H_ocean;  % [Pa]  Pressure at the ocean floor
P_frac = 2;             % [MPa] Pressure required to break rocks
P_floor = P_ocean*1e-6 + P_frac;    % [MPa] Total pressure to frac


% Piston/Wave Parameters
xi = 1;                 % [m]   WEC Amplitude
omega = 0.7;            % [rad/s]   Wave Frequency
v_amp = xi*omega;       % [m/s] Velocity Amplitude
piston_area = 0.41;     % [m^2] Area of Piston

% Pipe Parameters
H = H_ocean + H_ground; % [m]   Total length of pipe
d_pipe = 0.26;          % [m]   Dimater of pipe
A_pipe = d_pipe^2*pi/4; % [m^2] Area of pipe
N = 40;                 % [-]   Number of pipe segments

out = sim('CASHEW.slx');
plot(out.simout.Force_on_WEC)