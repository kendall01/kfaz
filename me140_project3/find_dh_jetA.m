function [ dh ] = find_dh_jetA( T1,T2,phi )
%Returns kJ

f_temp = @(x,y) sp_heats_jetA(x,y); %Note: Cp is the first element in the return of sp_heats and this returns just cp as a 1x1 double
dh = zeros(length(phi),1);
for i = 1:length(phi)
    dh(i) = integral( (@(x)f_temp(x,phi(i))) ,T1,T2);
end           
end

