function [ U ] = findU( T0,P0,mdot,A )
    R = 287;
    [~,~,k] = sp_heats(T0);
    
    % Find MFP via Eq. (1)
    MFP = (mdot./A).*sqrt(R.*T0)./P0;

    % Determine Ma that satisfies Eq. (2)
    M = sym('m',[length(MFP),1]);
    assume(M, 'real');
    assumeAlso(M < 1); %Select the subsonic solution. Set M > 1 if want supersonic solution.
    Ma = vpasolve( MFP == M.*sqrt(k).*(1+((k-1)./2).*M.^2).^(-(k+1)./(2.*(k-1))), M );
    Ma = struct(Ma);
    Ma = struct2array(Ma)'; %Ma is returned as a symbolic variable which doesn't work in polyval in sp_heats
    Ma = double(Ma);
    
    T = T0./(1 + (Ma.^2).*((k-1)./2));
    c = sqrt(k.*R.*T);
    U = Ma.*c;
end

