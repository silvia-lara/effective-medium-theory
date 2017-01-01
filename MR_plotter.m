%% MR Calculator
function[MR] = MR_plotter(eta, mu, B)

% eta = n/nrms; mean number of impurities over root mean square of the
% distribution of impurities

% Rho at a certain B field
[a,b]= sigmaemt(eta, mu, B);
rho_xx = a./(a.^2+b.^2);

%Rho at zero B field
[a0,b0]=sigmaemt(eta,mu, 0);
rho_0 = a0/(a0^2+b0^2);

%Magnetoresistance
MR=(rho_xx-rho_0)/rho_0;
end


