function [cp_m, cv_m, gamma_m] = sp_heats_mix(T,AF)
% Note: 
% (i) For this entire code refer to each species in the
% --- produces by these numbers:
% --- (1) CO2, (2) H2O, (3) N2, (4) O2
% (ii) Joules
% (iii) code inputs matrices

R = 287;                        % [J/kg*K]
N_to_O = 79/21;                 % Engineering Air Molar Mass Ratio of Nitrogen to Oxygen
AFs = 14.43;                    % stoichiometric Air-Fuel-Ratio, found from balancing equation (LEC 5, slides 3-4)
phi = AFs./AF;                  % Equivalence Ratio                                                  

a = [22.26 32.24 28.9 25.48];
b = [5.981*10^-2 0.1923*10^-2 -0.1571*10^-2 1.52*10^-2];
c = [-3.501*10^-5 1.055*10^-5 0.8081*10^-5 -0.7155*10^-5];
d = [7.469*10^-9 -3.595*10^-9 -2.873*10^-9 1.312*10^-9];

% Molar Mass of Each Species
MC = 12;
MH = 1;
MO = 16;
MN = 14;
M = [MC+2*MO 2*MH+MO 2*MN 2*MO];    % mores sig figs?
Mtotal = sum(M);

% Moles of Each Species (LEC 5,slide 8 [ASSUME: fuel-lean] )
N = [ phi, 2*phi, linspace(2*N_to_O, 2*N_to_O, length(phi))', 2*(1-phi) ]; 
Ntotal = sum(N);

% Mole Fraction of Each Species

for i = 1: length(N(1,:))
    y(:, i) = N(:,i)./Ntotal(i);
    mf(:,i) = ( y(:,i)*M(i) )/Mtotal;
end



P = zeros(4,4);
P(1,:) = [a(1) b(1) c(1) d(1)];
P(2,:) = [a(2) b(2) c(2) d(2)];
P(3,:) = [a(3) b(3) c(3) d(3)];
P(4,:) = [a(4) b(4) c(4) d(4)];

cp_m = 0;
cv_m = 0;
gamma_m = 0;
for j = 1:length(P)
% Find Cp for each Species
throw == error
cp = polyval(P,T);
cp = cp ./ M(j); %convert from KJ/kmol-K to J/kg-K
cv = cp - (R); %R converted to J/kg-K

% Find

cp_m = cp_m + mf(i)*cp(i);
cv_m = cv_m + mf(i)*cv(i);
gamma_m = cp_m./cv_m;
end




end