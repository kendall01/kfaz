function [ h_jetA, Tp ] = findAdiabaticFlameTemp( h_co2, h_h20, h_n2, h_o2, AF, LHV, Mfuel, TR )
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
error = 1E-4;
diff = 0.01;
speedFactor = 10;
N_to_O = 79/21;                 % Engineering Air Molar Mass Ratio of Nitrogen to Oxygen


% Find Formation Enthalpy (ASSUME: phi = 1, LEC 5,SLIDE 24) 
% ---------------------------------------------------------
dh_jetA = find_dh_mix(T0, TR, AF); %not correct, needs to refer to a new specific heats function that goes through the molar masses of the species in jet fuel
Q = LHV*Mfuel;
hf_jetA = -dh_jetA + 12.3*h_co2 + 11.1*h_h2o + Q; % enthaply of formation, jetA

% Find Adiabatic Flame Temperature (ASSUME: Q = 0, LEC 5,SLIDE 8)
% ---------------------------------------------------------------
% 5 SPECIES: (1)jet A, (2)co2, (3)h20, (4)n2, (5)o2 

%Function definitions for integral helpers
% ------------------------------
% TODO: WRITE THESE 4 FUNCTIONS!
% ------------------------------


% Enthalpies of Formation of the Products
hf = [h_co2, h_h20, h_n2, h_o2 ];

phi = linspace(0.05, 0.65);
% Find Molar Flow Rates
N = [phi, 2*phi, linspace(2*N_to_O, 2*N_to_O, length(phi))', 2*(1-phi) ]; 
sum = N * hf'; %intentional matrix multiplication multiplies and sums as desired.

% Solve for Tp (ASSUME: Q = 0, adiabatic)
% LHS: hf_jetA + dh_jetA - sum( N(i)*hf(i) )
% RHS: integral[ sum( Cp(i)(T')dT') ] from T0 to Tp + Q
LHS = hf_jetA + dh_jetA - sum;

for p = 0:length(phi)


    % Initializations
    T1 = T0;
    T2 = T1 + .01;
    
    diff = RHS(T1,T2) - LHS; %not correct, needs n_dot
    iterations = 0;
    
    % ----------------------------------------------------------------
    % TODO: ADD FOR LOOP?, NEEDS TO GO THROUGH MULTIPLE TEMPERATURES!!
    % -----------------------------------------------------------------
    % ( modeled after applyIsentropicTempVar from Project #1 )
    % for
    while(abs(diff) > error)
        T2 = T2 - diff./speedFactor;
        RHS = find_dh_mix(T1,T2,phi(p));
        diff = RHS - LHS;
        iterations = iterations + 1;
    end
    iterations
    % end
end
    
end



