% PsatW.m 
% 4-22-16 - Created Jon Renslo
function [psat] = PsatW(T)
    %psat in pascals
    psat = exp(-1.2914*10^8./T.^3 +8.2048*10^5./T.^2 -6522.8./T +25.5887);
end
