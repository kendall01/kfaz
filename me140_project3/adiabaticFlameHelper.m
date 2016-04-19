function [RHS] = adiabaticFlameHelper(T1, T2)
f_temp_co2 = @(x) sp_heats_co2(x); 
f_temp_h2o = @(x) sp_heats_h2o(x);
f_temp_n2  = @(x) sp_heats_n2(x);
f_temp_o2  = @(x) sp_heats_o2(x);

cp_co2 = @(a,b) integral(f_temp_co2,a,b);
cp_h2o = @(a,b) integral(f_temp_h2o,a,b);
cp_n2  = @(a,b) integral(f_temp_n2,a,b);
cp_o2 = @(a,b) integral(f_temp_o2,a,b);

RHS = @(a,b) cp_co2(a,b) + cp_h2o(a,b) + cp_n2(a,b) + cp_o2(a,b);