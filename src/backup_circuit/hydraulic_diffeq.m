

[tout,yout] = ode45(dynamics,tspan,y0);


[dState_dt, outputs] = function dynamics(t,state)

rho = state(1); % density of CO2
V = state(2); % volume of CO2 in piston
Vdot = state(3); % time derivative of volume of CO2 in piston



end
