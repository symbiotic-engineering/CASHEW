function [T_CO2, T_s_in, T_s_out, T_insu] = temp_model(P_vs_depth,rho_vs_depth,...
                    outer_radius_pipe,inner_radius_pipe,outer_radius_insu,...
                    thickness_pipe,thickness_insulation,T_supercritical,T_oc,...
                    k_pipe,k_insu,h_in,h_out,P_heat,rho_pipe,c_CO2,c_pipe,...
                    h_under_seafloor,deltaZ,g,N)

deltaP = [0 diff(P_vs_depth)];
rho_CO2 = rho_vs_depth;

dA_lat_pipe_outer = 2 * pi * outer_radius_pipe * deltaZ; % lateral area of a sliver of outside of pipe with length deltaZ, m^2
dA_lat_pipe_inner = 2 * pi * inner_radius_pipe * deltaZ; % lateral area of a sliver of inside  of pipe with length deltaZ, m^2
dA_lat_insu_outer = 2 * pi * outer_radius_insu * deltaZ; % lateral area of a sliver of outside of insulation with length deltaZ, m^2

A_cross_pipe = pi * (outer_radius_pipe^2 - inner_radius_pipe); % cross sectional area of pipe, m^2
A_cross_CO2 = pi * inner_radius_pipe^2;                        % cross sectional area of CO2,  m^2

heat_capacity_CO2  = A_cross_CO2  * deltaZ * rho_CO2   * c_CO2;  % mc for a sliver of CO2  with length deltaZ, J/C
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
T_CO2   = zeros(1,N); % temperature of CO2

%initialize TCO2 at the surface of the ocean
T_CO2(1) = T_supercritical;

for i = 1:N-1
    heat_xfer_from_CO2_to_ocean = 1/total_thermal_resistance * (T_CO2(i) - T_oc);
    net_heat_xfer_into_CO2  = -heat_xfer_from_CO2_to_ocean + g * deltaZ + deltaP(i) / rho_CO2(i) + P_heat;

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

end