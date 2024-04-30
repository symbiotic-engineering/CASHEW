piston_area = 0.26;     %   [m^2]   Area of Piston
rhoTPData = readmatrix("../../data/CO2TempData.csv");
temp = [30;40];
pressure = flip([30;20;10;9;8;7;6;5;4;3;2;1;0.007]);
rhoAt30 = rhoTPData(:,2); 
rhoAt40 = rhoTPData(:,5); 
rho = [rhoAt30 rhoAt40]';

