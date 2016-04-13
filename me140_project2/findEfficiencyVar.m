function [ n ] = findEfficiencyVar( T1,T2,T2s )

error = 1E-4;
dT = 0.01;
speed_ratio = 1000;
f_temp = @(x) sp_heats(x); %Note: Cp is the first element in the return of sp_heats and this returns just cp as a 1x1 double
n = integral( f_temp,T1,T2s )/integral( f_temp,T1,T2 );
             
end


h