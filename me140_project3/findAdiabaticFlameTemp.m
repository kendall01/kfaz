function [ Tp ] = findAdiabaticFlameTemp( h_co2 h_h20 )
%THINGS TO DO:
% 1) put actual integral into line 20 (Check Red Rectangle in pdf)
% 2) add in section for finding adiabatic flame temperature

% INPUTS:
% --------
%Find h_co2 and h_h2o from Tables in units of KJ/mol
LHV = 42800; %kJ/kg
Mfuel = 170 kg/mol

% OUTPUTS:
% ---------
% Tp adiabatic flame temperature

% ASSUME: 
% --------
%

% Find Formation Enthalpy: 
% ---------

h_jetA = 12.3*h_co2 + 11.1*h_h2o + (LHV*M_fuel) - %integral of Cp_JetA*dT from Tr to T_ref


% Find Adiabatic Flame Temperature
% ---------
phi = linspace(0.05, 0.65);
for i = 0:length(phi)
    N_jetA(i) = phi(i);
    N_co2(i) = phi(i);
    N_h2o(i) = 2*phi(i);
    N_N2(i) = 2*3.76;
    N_o2(i) = 2*(1-phi(i));
    
    %Solve for Tp
end




end

