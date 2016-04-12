function [ T,T0,Ma,U,k ] = findMaTemps( mdot,A,Tm,P0,RF )
% ASSUME: Constant Cp,Cv - because we are operating at low Mach numbers, 
% we can assume the difference between assuming constant gamma and varrying
% gamma is small. 
% FIND: Ma, T & T0, U (velocity), k (gamma)
% NOTE: Approach from ME140 Lecture 3 page 20



error = 0.01;
R = 287;

% (a) Assume: T0guess = Tm 
% Initializations
T0guess = Tm;               
k = sp_heats(T0guess);     
i = 1;

while i==1;
    % (b) Find MFP via Eq. (1)
    MFP = (mdot/A)*sqrt(R*T0guess)/P0;

    % (c) Determine Ma that satisfies Eq. (2)
    syms M 
    Ma = solve( MFP == M*sqrt(k)*(1+((k-1)/2)*M^2)^(-(k+1)/(2*(k-1))), M );

    % (d) Find T via Eq (3)
    T = Tm/(1+RF*((k-1)/2)*Ma^2);

    % (e) Find T0 via Eq. (4), If T0 = T0(initial guess) then T0 is correct.
    T0 = T*(1+((k-1)/2)*Ma^2);

    if abs(T0-T0guess)<error
        i = 0;
    else
        % --- Else, let T0guess = T0 and repeat from (b)
        T0guess = T0;
    end 
    k = sp_heats(T); 
end

U = Ma*sqrt(gamma*R*T0);
end

