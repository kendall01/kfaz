%hEng.m , wrapper for energyF
% 4-22-16 created by Jon Renslo
function H = hEng(T,spec,mol)
    
    if(nargin == 3)
        mol = 1;
    end
    a = energyF(T,1e5,spec,mol); %pressure does not affect H
    H = a.H;
end