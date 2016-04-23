%gEng.m , wrapper for energyF
% 4-22-16 created by Jon Renslo
function G = gEng(T,P,spec,mol)
    
    if(nargin == 3)
        mol = 1;
    end
    a = energyF(T,P,spec,mol);
    G = a.G;
end