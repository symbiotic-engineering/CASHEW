piston_area = 0.26;     %   [m^2]   Area of Piston
rhoTPData = readmatrix("../../data/CO2TempData.csv");
temp = [30;40]; % degC
pressure = flip([30;20;10;9;8;7;6;5;4;3;2;1;0.007]); % MPa
rhoAt30 = rhoTPData(:,2); 
rhoAt40 = rhoTPData(:,5); 
rho = [rhoAt30 rhoAt40]; % kg/m^3

dP = [diff(pressure)]; 
dP = [dP; dP(length(dP))]; % append final value

dRho = [diff(rho)];
dRho = [dRho ; dRho(length(dRho),:)]; % append final row

bulkModulus = rho .* ([dP dP]./dRho);
bulkModulus = bulkModulus'; %  adjust for input to piston
rho = rho'; % adjust for input to piston

cp = ones(size(rho))*0.846;

dRhoP = diff(rho);
dRhoP = [dRhoP;dRhoP];
for i=1:3
    for j=1:2
        dRhoP(j,i) = -0.1;
    end
end

dT = ones(size(rho))*diff(temp);
alpha = abs(-(1./rho) .* (dRhoP./dT));

tempMatrix = repmat(temp,1,length(rho));
cv= cp - (tempMatrix .* bulkModulus .* (alpha).^2)./rho;

