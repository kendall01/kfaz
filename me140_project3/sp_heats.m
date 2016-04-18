function [cp,cv,gamma] = sp_heats(T)
%Works for matrices
% Returns Joules
a = 28.11;
b = 0.1967E-2;
c = 0.4802E-5;
d = -1.966E-9;
molar_mass_air = .02897;
R = 287;                        % [J/kg*K]

P = [d,c,b,a];

cp = polyval(P,T);
cp = cp ./ molar_mass_air; %convert from KJ/kmol-K to J/kg-K
cv = cp - (R); %R converted to J/kg-K
gamma = cp./cv; 

end