function [ dh ] = find_dh( T1,T2 )
%Returns kJ
error = 1E-4;
dT = 0.01;
speed_ratio = 1000;
f_temp = @(x) sp_heats(x); %Note: Cp is the first element in the return of sp_heats and this returns just cp as a 1x1 double
dh = zeros(length(T1),1);
for i = 1:length(T1)
    dh(i) = integral( f_temp,T1(i),T2(i) );
end           
end

