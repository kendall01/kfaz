function [ T2 ] = applyIsentropicTempVar_mix( T1, P2_by_P1, phi)
% DESCRIPTION: 
% INPUTS: T1, inlet temperature
% ------- P2_by_P1, ratio between inlet and exit pressures
% ------- k, Cp/Cv ratio of specific heats
% OUTPUT: T2, exit temperature
% ASSUME: isentropic
error = 1E-1;
R = 287;                        % [J/kg*K]

%Function definitions for integral helpers, etc.
f_temp = @(x) sp_heats_mix(x,phi)./x; %Note: Cp is the first element in the return of sp_heats and this returns just cp as a 1x1 double
RHS = @(a,b) integral(f_temp,a,b);

%Array initializations
T2 = zeros(1,length(T1));

parfor i = 1:length(T1)
    LHS = R*log(P2_by_P1(i));
    T2(i) = T1(i)+.01;
    diff = RHS(T1(i),T2(i)) - LHS;
    iterations = 0;
    while(abs(diff) > error)
        T2(i) = T2(i) - diff./1E0;
        diff = RHS(T1(i),T2(i)) - LHS;
        iterations = iterations + 1;
    end
    iterations;
end

end


