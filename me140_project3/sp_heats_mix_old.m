function [cp_m, cv_m, gamma_m] = sp_heats_mix_old(T,AF)
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
Ntotal = sum(N,2); %Sums along vertical dimension



% Mole Fraction of Each Species

for i = 1: size(N,2) %number of columns in N
    y(:, i) = N(:,i)./Ntotal;
    mf(:,i) = ( y(:,i)*M(i) )/Mtotal;
end

P = zeros(4,4);
P(1,:) = [a(1) b(1) c(1) d(1)];
P(2,:) = [a(2) b(2) c(2) d(2)];
P(3,:) = [a(3) b(3) c(3) d(3)];
P(4,:) = [a(4) b(4) c(4) d(4)];

cp_m = zeros(size(T));
cv_m = zeros(size(T));
gamma_m = zeros(size(T));

for col = 1:size(T,2)

    for row = 1:size(T,1) %Size returns number of rows in T
        for spc = 1:length(P)
        % Find Cp for each Species
        cp = polyval(P(spc,:),T(row, col));
        cp = cp ./ M(spc); %convert from KJ/kmol-K to J/kg-K
        cv = cp - (R); %R converted to J/kg-K

        % Find Cp for Mixture
        cp_m(row, col) = cp_m(row, col) + mf(row,spc)*cp;
        cv_m(row, col) = cv_m(row, col) + mf(row,spc)*cv;
        gamma_m(row, col) = cp_m(row, col)./cv_m(row, col);
        end

    end
    
end




end