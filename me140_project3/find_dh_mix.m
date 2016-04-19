function [ dh ] = find_dh_mix( T1,T2,AF )
%Returns kJ

f_temp = @(x,y) sp_heats_mix(x,y); %Note: Cp is the first element in the return of sp_heats and this returns just cp as a 1x1 double
dh = zeros(size(T1,1),1);
for i = 1:size(T1,1)
    dh(i) = integral( (@(x)f_temp(x,AF(i))) ,T1(i),T2(i));
end           
end

