function [ w x y z ] = findProductMix( mair mfuel Mfuel )
% INPUTS:
% --------
% mair, mass of air
% mfuel, mass of fuel
% Mfuel, molar mass of fuel

% OUTPUTS:
% ---------
% w,x,y,x = mols of CO2, H20, N2, O2

% ASSUME: 
% --------
% (i) complete combustion 
% (ii) variable Cp,Cv

% COMBUSTION REACTION: 
% -------------------
% N_fuel*(fuel) + N_air*(02 + 3.76N2) --> w*CO2 +x*H20 + y*N2 +z*O2

% CONSTANTS
Mair = 28.97; %[g/mol]

Nair = mair./Mair;
Nfuel = mfuel./Mfuel;

a= 
b= 
c= 

ac = [a*Nfuel e*Nair 0 0 0];
ah = [b*Nfuel f*Nair 0 0 0];
ao = [c*Nfuel g*Nair 0 0 0];
an = [d*Nfuel h*Nair 0 0 0];






end

