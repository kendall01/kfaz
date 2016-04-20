function [cp_m, cv_m, gamma_m, cp_co2, cp_h2o, cp_n2, cp_o2] = sp_heats_jetA(T,phi)
% Note: 
% (i) For this entire code refer to each species in the
% --- produces by these numbers:
% --- (1) CO2, (2) H2O, (3) N2, (4) O2
% (ii) Joules
% (iii) code inputs matrices of temperature. should work for multiple rows
% and multiple columns simultaneously. definitely works for either
% independently.

R = 287;                        % [J/kg*K]
N_to_O = 79/21;                 % Engineering Air Molar Mass Ratio of Nitrogen to Oxygen

a = [22.26 32.24 28.9 25.48]; %[kJ/kmol-K]
b = [5.981*10^-2 0.1923*10^-2 -0.1571*10^-2 1.52*10^-2];
c = [-3.501*10^-5 1.055*10^-5 0.8081*10^-5 -0.7155*10^-5];
d = [7.469*10^-9 -3.595*10^-9 -2.873*10^-9 1.312*10^-9];

% Molar Mass of Each Species
MC = 12; %[g/mol]
MH = 1;
MO = 16;
MN = 14;
M = [MC+2*MO 2*MH+MO 2*MN 2*MO];    % mores sig figs?
Mtotal = sum(M); %[g/mol]

% Moles of Each Species (LEC 5,slide 8 [ASSUME: fuel-lean] )
N = [ phi, 2*phi, linspace(2*N_to_O, 2*N_to_O, length(phi))', 2*(1-phi) ]; 
Ntotal = sum(N,2); %Sums along vertical dimension


% Mole Fraction of Each Species

for i = 1: size(N,2) %number of columns in N
    y(:, i) = N(:,i)./Ntotal;
    mf(:,i) = ( y(:,i)*M(i) )/Mtotal;
end

P = zeros(4,4);
P(s.CO2,:) = [a(1) b(1) c(1) d(1)];
P(s.H2O,:) = [a(2) b(2) c(2) d(2)];
P(s.N2,:) = [a(3) b(3) c(3) d(3)];
P(s.O2,:) = [a(4) b(4) c(4) d(4)];

cp_m = zeros(size(T));
cv_m = zeros(size(T));
gamma_m = zeros(size(T));
cp_co2 = zeros(size(T));
cp_h2o = zeros(size(T));
cp_n2 = zeros(size(T));
cp_o2 = zeros(size(T));

for col = 1:size(T,2)

    for row = 1:size(T,1) %Size returns number of rows in T
        for spc = 1:length(P)
        % Find Cp for each Species
        cp = polyval(P(spc,:),T(row, col));
        cp = cp ./ M(spc); %convert from KJ/kmol-K to J/kg-K
        cv = cp - (R); %R converted to J/kg-K

        switch spc
            case s.CO2
                cp_co2(row,col) = cp;
            case s.H2O
                cp_h2o(row,col) = cp;
            case s.N2
                cp_n2(row,col) = cp;
            case s.O2
                cp_o2(row,col) = cp;
            otherwise
                throw == error;
        end
            
        
        % Find Cp for Mixture
        cp_m(row, col) = cp_m(row, col) + mf(row,spc)*cp;
        cv_m(row, col) = cv_m(row, col) + mf(row,spc)*cv;
        gamma_m(row, col) = cp_m(row, col)./cv_m(row, col);
        end

    end
    
end




end