function [ T ] = recoveryFactor( Tm,Ma,RF )
k = sp_heats(Tm);
T  = Tm/(1+RF*((k-1)/2)*Ma^2); 


end

