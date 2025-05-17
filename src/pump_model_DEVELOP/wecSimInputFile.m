%% Simulation Data
simu = simulationClass();                % Initialize simulationClass
simu.simMechanicsFile = 'CASHEW.slx';    % Simulink Model File
simu.startTime = 0;                      % Simulation Start Time [s] 
simu.rampTime = 0;                       % Wave Ramp Time [s] 
simu.endTime= 100;                       % Simulation End Time [s]
simu.solver = 'daessc';                  % Solver; 'ode45', 'ode23tb','daessc', etc. daessc seems to work best
simu.mode = 'normal';                    % 'accelerator', 'normal' 
simu.dt = 0.01;                         % Simulation time-step [s] 0.001 may give slightly better results.

%% Wave Information  
% Regular Waves 
waves = waveClass('regular');   % Initialize waveClass
waves.height = 2;               % Wave Height [m]
waves.period = 8.97597901;      % Wave Period [s]

%% Body Data
% Float
body(1) = bodyClass('rm3.h5');                           % Initialize bodyClass for Float
body(1).geometryFile = 'float.stl';                      % Geometry File
body(1).mass = 727010;                            % Mass [kg] - actual mass of 727010 kg; can replace 'equilibrium'
body(1).inertia = [20907301 21306090.66 37085481.11];    % Moment of Inertia [kg*m^2] 

% Spar/Plate
body(2) = bodyClass('rm3.h5');                             % Initialize bodyClass for Spar
body(2).geometryFile = 'plate.stl';                        % Geometry File
body(2).mass = 878300;                              % Mass [kg] - actual mass of 878300 kg; can replace 'equilibrium'
body(2).inertia = [94419614.57 94407091.24 28542224.82];   % Moment of Inertia [kg*m^2]

%% PTO and Constraint Parameters
% Floating (3DOF) Joint
constraint(1) = constraintClass('Constraint1');   % Initialize constraintClass for Constraint1
constraint(1).location = [0 0 0];                 % Constraint Location [m]

% Translational PTO
pto(1) = ptoClass('PTO1');   % Initialize ptoClass for PTO1
pto(1).stiffness = 0;        % PTO Stiffness [N/m]
pto(1).damping = 0;          % PTO Damping [N/(m/s)]
pto(1).location = [0 0 0];   % PTO Location [m]