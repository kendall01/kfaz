function [ T,T0,Ma,U,k ] = findMaTemps( mdot,A,Tm,P0,RF )
% ASSUME: Constant Cp,Cv - because we are operating at low Mach numbers, 
% we can assume the difference between assuming constant gamma and varrying
% gamma is small. 
% FIND: Ma, T & T0, U (velocity), k (gamma)
% NOTE: Approach from ME140 Lecture 3 page 20
% INPUT: Accepts scalar or vector inputs. Vector inputs for a two state
% system consisting of state a and state b should be entered as column
% vectors as illustrated: findMaTemps([mdot_a; mdot_b], [A_a; A_b], [Tm_a;
% Tm_b], [P0_a; P0_b], [RF_a; RF_b])

error = 0.01;
R = 287;

% (a) Assume: T0guess = Tm 
% Initializations
T0guess = Tm;
T0 = T0guess + 2 * error;
[~,~,k] = sp_heats(T0guess);     

while abs(T0-T0guess)>error;
    % (b) Find MFP via Eq. (1)
    T0guess = T0;

    MFP = (mdot./A).*sqrt(R.*T0guess)./P0;

    % (c) Determine Ma that satisfies Eq. (2)
    M = sym('m',[length(MFP),1]);
    assume(M, 'real');
    assumeAlso(M < 1); %Select the subsonic solution. Set M > 1 if want supersonic solution.
    Ma = vpasolve( MFP == M.*sqrt(k).*(1+((k-1)./2).*M.^2).^(-(k+1)./(2.*(k-1))), M );
    Ma = struct(Ma);
    Ma = struct2array(Ma)'; %Ma is returned as a symbolic variable which doesn't work in polyval in sp_heats
    Ma = double(Ma);
    
    % (d) Find T via Eq (3)
    T = Tm./(1+RF.*((k-1)./2).*Ma.^2);

    % (e) Find T0 via Eq. (4), If T0 = T0(initial guess) then T0 is correct.
    T0 = T.*(1+((k-1)./2).*Ma.^2);

    % update
    [~,~,k] = sp_heats(T); 
end

U = Ma.*sqrt(k.*R.*T0);
end

