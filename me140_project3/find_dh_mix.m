function [ dh ] = find_dh_mix( T1,T2,AF )
%Returns kJ
AFs = 14.43;                    % stoichiometric Air-Fuel-Ratio, found from balancing equation (LEC 5, slides 3-4)

f_temp = @(x,y) sp_heats_mix(x,y); %Note: Cp is the first element in the return of sp_heats and this returns just cp as a 1x1 double
dh = zeros(size(T1,1),1);
for i = 1:size(T1,1)
    dh(i) = integral( (@(x)f_temp(x,AFs./AF(i))) ,T1(i),T2(i));
end           
end

