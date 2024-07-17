% example function call: get_P_surface(1,2700,2900,1000*9.8*2700,9.8,1,1e-17)

function [P_surface_required,P_bottom, P_vs_depth, rho_vs_depth] = get_P_surface(mdot,h_seafloor,h_under_seafloor,P_bottom_required,g,D,k)

    % iterate until P_surface solves the nonlinear equation
    P_surface_guess = 8e6; % 8 MPa
    P_bottom_error_fcn = @(P_surface_guess) (pressure_vs_depth_fcn(mdot,P_surface_guess,h_seafloor,h_under_seafloor,g,D,k) - P_bottom_required);
    P_surface_required = fzero(P_bottom_error_fcn, P_surface_guess);

    % rerun result with converged P_surface to get depth profile
    [P_bottom, P_vs_depth, rho_vs_depth] = pressure_vs_depth_fcn(mdot,P_surface_required,h_seafloor,h_under_seafloor,g,D,k);

end

function [P_bottom, P_vs_depth, rho_vs_depth] = pressure_vs_depth_fcn(mdot,P_surface,h_seafloor,h_under_seafloor,g,D,k)
    rhoCO2ref = [0 10 30 50 85 110 145 235 500 620 700 790 860 905 935]; % kg/m^3
    pCO2ref   = [0 1  2  3  4  5   6   7   8   9   10  15  20  25  30 ]*1e6; % Pa
    % this data assumes temperature of approx 35 degree C

    P_vs_depth = zeros(1,h_under_seafloor);
    P_vs_depth(1) = P_surface;
    rho_vs_depth = zeros(1,h_under_seafloor);
    
    %% iteration to calculate CO2 pressure at depth
    for i=2:h_under_seafloor
        rho_i = interp1(pCO2ref, rhoCO2ref, P_vs_depth(i-1),'linear','extrap'); % density of CO2, kg/m^3
        mu_i = 3e-6; % dynamic viscosity of CO2 [Pa s]

        % note: delta_h_i=1 by convention
        deltaP_hydro_i = g * rho_i; % Pa, iterative hydrostatic pressure

        Q = mdot / rho_i; % volume flow rate
        A = pi/4 * D^2;   % area
        v = Q/A;          % velocity

        if i < h_seafloor % pipe loss
            f = 0.015; % turbulent friction factor, taken for roughness ~5e-4 at high Re on Moody chart.
            deltaP_loss_i = 1/2 * rho_i * v^2 * f/D;
        else % injection loss
            fracking = true;
            if fracking
                P_fracking = 2e6 / 200; % hardcoded for now, corresponding to 2MPa over 200m
                deltaP_loss_i = P_fracking;
            else % linear Darcy loss based on rock permeability
                deltaP_loss_i = mu_i * v / k;
            end
        end
        deltaP_i = deltaP_hydro_i - deltaP_loss_i;

        rho_vs_depth(i) = rho_i;
        P_vs_depth(i) = P_vs_depth(i-1) + deltaP_i;
    end
    P_bottom = P_vs_depth(end);

    rho_vs_depth(1) = rho_vs_depth(2);

end
