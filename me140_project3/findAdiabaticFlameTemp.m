function [ h_jetA, Tp ] = findAdiabaticFlameTemp( h_co2, h_h20, h_n2, h_o2, AF, LHV, Mfuel )
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
TR = 300;                       % Temperature of Reactants, [K]
error = 1E-4;
diff = 0.01;
speedFactor = 10;

% Find Formation Enthalpy (ASSUME: phi = 1, LEC 5,SLIDE 24) 
% ---------------------------------------------------------
dh_jetA = find_dh_mix(T0, TR, AF);
Q = LHV*Mfuel;
hf_jetA = -dh_jetA + 12.3*h_co2 + 11.1*h_h2o + Q; % enthaply of formation, jetA

% Find Adiabatic Flame Temperature (ASSUME: Q = 0, LEC 5,SLIDE 8)
% ---------------------------------------------------------------
% 5 SPECIES: (1)jet A, (2)co2, (3)h20, (4)n2, (5)o2 

%Function definitions for integral helpers
% ------------------------------
% TODO: WRITE THESE 4 FUNCTIONS!
% ------------------------------
f_temp_co2 = @(x) sp_heats_co2(x); 
f_temp_h2o = @(x) sp_heats_h2o(x);
f_temp_n2  = @(x) sp_heats_n2(x);
f_temp_o2  = @(x) sp_heats_o2(x);

cp_co2 = @(a,b) integral(f_temp_co2,a,b);
cp_h2o = @(a,b) integral(f_temp_h2o,a,b);
cp_n2  = @(a,b) integral(f_temp_n2,a,b);
cp_o2 = @(a,b) integral(f_temp_o2,a,b);

RHS = @(a,b) cp_co2(a,b) + cp_h2o(a,b) + cp_n2(a,b) + cp_o2(a,b);

% Enthalpies of Formation
hf = [ hf_jetA, h_co2, h_h20, h_n2, h_o2 ];

numSpecies = 5;                     
phi = linspace(0.05, 0.65);
N = zeros(length(phi),numSpecies);
for p = 0:length(phi)
    % Find Molar Flow Rates
    N(p,1) = phi(p);
    N(p,2) = phi(p);
    N(p,3) = 2*phi(p);
    N(p,4) = 2*3.76;
    N(p,5) = 2*(1-phi(p));
    
    for s = 1:numSpecies
        sum = sum + N(p,s)*hf(s);
    end
    
    % Solve for Tp (ASSUME: Q = 0, adiabatic)
    % LHS: hf_jetA + dh_jetA - sum( N(i)*hf(i) )
    % RHS: integral[ sum( Cp(i)(T')dT') ] from T0 to Tp + Q
    LHS = hf_jetA + dh_jetA - sum;
    
    % Initializations
    T1 = T0;
    T2 = T1 + .01;
    
    RHS = RHS(T1,T2);
    diff = RHS - LHS;
    iterations = 0;
    
    % ----------------------------------------------------------------
    % TODO: ADD FOR LOOP?, NEEDS TO GO THROUGH MULTIPLE TEMPERATURES!!
    % -----------------------------------------------------------------
    % ( modeled after applyIsentropicTempVar from Project #1 )
    % for
    while(abs(diff) > error)
        T2 = T2 - diff./speedFactor;
        diff = RHS(T1,T2) - LHS;
        iterations = iterations + 1;
    end
    iterations
    % end
end
    
end



