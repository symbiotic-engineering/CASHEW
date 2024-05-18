piston_area = 0.26;     %   [m^2]   Area of Piston
initialRead = readmatrix("../../data/CO2TempData.csv");
[r,c] = size(initialRead);
temp = initialRead(1,2:c); % degC
pressure = initialRead(2:r,1); % MPa
rho = initialRead(2:r,2:c); % kg/m^3

dP = [diff(pressure)]; 
dP = [dP; dP(length(dP))]; % append final value
dP = repmat(dP,1,c-1);

dRho = [diff(rho)];
dRho = [dRho ; dRho(length(dRho),:)]; % append final row

bulkModulus = rho .* (dP./dRho);
bulkModulus = bulkModulus'; %  adjust for input to piston
rho = rho'; % adjust for input to piston

cp = ones(size(rho))*1.846; %%%%%%%%%% changed to eliminate issue

dRhoP = diff(rho);
dRhoP = [dRhoP; dRhoP(c-2,:)];
for i=1:r-1
    for j=1:c-1
        if dRhoP(j,i) == 0
            dRhoP(j,i) = -0.1;
        end
    end
end

dT = diff(temp);
dT = [dT dT(length(dT))];
dT = repmat(dT,r-1,1)';

alpha = abs(-(1./rho) .* (dRhoP./dT));

% tempMatrix = repmat(temp,1,length(rho))+273.15;
% cv = cp - (tempMatrix .* bulkModulus .* (alpha).^2)./rho;

