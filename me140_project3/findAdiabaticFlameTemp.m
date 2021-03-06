function [ hf_jetA, T2f ] = findAdiabaticFlameTemp( h_co2, h_h2o, h_n2, h_o2, LHV, Mfuel, TR, phi )
% INPUTS:
% --------
% h_co2 and h_h2o (from Tables in units of KJ/mol)

% OUTPUTS:
% ---------
% Tp adiabatic flame temperature
% h_jetA, enthalpy of formation of Jet A 

% ASSUME: 
% --------
% 

% CONSTANTS

CELSIUS_TO_KELVIN = 273.15;
T0 = 25 + CELSIUS_TO_KELVIN;    % Reference Temperature, [K] 
error = 1E2;
diff = 0.01;
speedFactor = 5E2;
N_to_O = 79/21;                 % Engineering Air Molar Mass Ratio of Nitrogen to Oxygen
G_to_KG = .001;

% Find Formation Enthalpy (ASSUME: phi = 1, LEC 5,SLIDE 24) 
% ---------------------------------------------------------
dh_jetA = 0; %assumption
Q = (LHV*G_to_KG)*Mfuel;
hf_jetA = -dh_jetA + 12.3*h_co2 + 11.1*h_h2o + Q; % enthaply of formation, jetA

% Find Adiabatic Flame Temperature (ASSUME: phi = linspace, Q = 0, LEC 5,SLIDE 8)
% ---------------------------------------------------------------
% 5 SPECIES: (1)jet A, (2)co2, (3)h20, (4)n2, (5)o2  - not really the case
% anymore as written

%Function definitions for integral helpers
% ------------------------------
% TODO: WRITE THESE 4 FUNCTIONS!
% ------------------------------


% Enthalpies of Formation of the Products
hf = [h_co2, h_h2o, h_n2, h_o2 ];

% Find Molar Flow Rates
N = [phi, 2*phi, linspace(2*N_to_O, 2*N_to_O, length(phi))', 2*(1-phi) ]; 
sum = N * hf'; %intentional matrix multiplication multiplies and sums as desired.

% Solve for Tp (ASSUME: Q = 0, adiabatic)
% LHS: hf_jetA + dh_jetA - sum( N(i)*hf(i) )
% RHS: integral[ sum( Cp(i)(T')dT') ] from T0 to Tp + Q
LHS = phi*hf_jetA + dh_jetA - [sum,sum];
T2f = zeros(length(phi),length(hf_jetA));
for l = 1:length(hf_jetA)
parfor p = 1:length(phi)
    % Initializations
    T1 = T0;
    T2 = T1 + .01;
    
    RHS = find_dh_mix(T1,T2,phi(p));
    diff = RHS - LHS(p,l); %not correct, needs n_dot
    iterations = 0;
    
    % modeled after applyIsentropicTempVar from Project #1
    while(abs(diff) > error)
        T2 = T2 - diff./speedFactor;
        RHS = find_dh_mix(T1,T2,phi(p));
        diff = RHS - LHS(p,l);
        iterations = iterations + 1;
    end
    iterations;
    T2f(p,l) = T2;
end
end
    
end



