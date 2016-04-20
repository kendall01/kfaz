% ME 140 PROJECT 3: SR30 TURBOJET ENGINE (COMBUSTION REPORT)
% FRANKIE WILLCOX, KENDALL FAGAN, ZACH NEVILLS, ANTOINE SCREVE
% -------------------------------------------------------------
% ASSUME: 
% (i) lean air conditions

clear all;
close all;

% UNIT CONVERSIONS
% ----------------
C_TO_K = 273.15;
LBF_TO_N = 4.44822;
IN2_TO_M2 = 6.4516*10^-4;
KPA_TO_PA = 10^3;
KG_TO_G = 10^3;
KJ_TO_J = 10^3;
R = 287.058;                                            % J/(kg*K)
plotfixing = 1;       % Boolean for whether or not to run plot fixer

% DATA:
% ------
rpm = me140_project2_data(1)';
mdot_fuel = me140_project2_data(13)';                   % mdot_fuel,[kg/s]
Fthrust_meas = me140_project2_data(14)' .* LBF_TO_N;    % Fthrust, [N]
Toil = me140_project2_data(7)' + C_TO_K; 

% Measured Temperatures
T2m = me140_project2_data(2)' + C_TO_K;                 % T2m, cross-flow
T3m = me140_project2_data(3)' + C_TO_K;                 % T3m, cross-flow 
T4m = me140_project2_data(4)' + C_TO_K;                 % T4m, cross-flow
T5m = me140_project2_data(5)' + C_TO_K;                 % T5m, axial-flow
T8m = me140_project2_data(6)' + C_TO_K;                 % T8m, cross-flow

% Measured Pressures
Patm = 101300;  
dP2 = me140_project2_data(8)'  .* KPA_TO_PA;            % dP2 = Pstag - Pstatic 
P03 = me140_project2_data(9)'  .* KPA_TO_PA + Patm ;    % P03 [absolute] [Pa]
P4 =  me140_project2_data(10)' .* KPA_TO_PA + Patm;    % P4 (static) [absolute] 
P05 = me140_project2_data(11)' .* KPA_TO_PA + Patm;    % P05 [absolute]   
P08 = me140_project2_data(12)' .* KPA_TO_PA + Patm;    % P08 [absolute]   

% GIVEN:
% ------
RFcross = linspace(.68, .68, length(rpm))';            % Reference Factor cross-section
RFaxial = linspace(.86, .86, length(rpm))';            % Reference Factor axial
hf_co2 =  -393.509 * KJ_TO_J;                          % [J/mol]
hf_h2o =  -241.818 * KJ_TO_J;
hf_n2  = 0         * KJ_TO_J;
hf_o2  = 0         * KJ_TO_J;

% Jet A (Kerosene) vs. Dodecane (C12H26)
HC_ratio = [1.8 2.17];
Mfuel = [170 170];
latent_heat_vap= [42800*KJ_TO_J 44560*KJ_TO_J ];        % latent heat of vaporization
LHV = [42800*KJ_TO_J, 44560*KJ_TO_J, 42900*KJ_TO_J ];   % vapor, J/kg Jet-A, Dodecane, TS-1 (source:Aviation Tech Review)
hf0 = [0 -1712];                                        % vapor, enthalpy of formation J/kg (replace 0 later)
AFs = 14.43;                                            % stoichiometric Air-Fuel-Ratio, found from balancing equation (LEC 5, slides 3-4)


% Areas [m^2]
A1 = linspace(27.3, 27.3, length(rpm))' .* IN2_TO_M2;  % Flow area at inlet of bellmouth
A2 = linspace(6.4, 6.4, length(rpm))'   .* IN2_TO_M2;  % Effective flow area at pitot-static probe
A3 = linspace(9.0, 9.0, length(rpm))'   .* IN2_TO_M2;  % Area of compressor exit
A4 = linspace(7.2, 7.2, length(rpm))'   .* IN2_TO_M2;  % Area before bend to turbine inlet
A5 = linspace(4.7, 4.7, length(rpm))'   .* IN2_TO_M2;  % Area of turbine outlet
A8 = linspace(3.87, 3.87, length(rpm))' .* IN2_TO_M2;  % Area of nozzle exit

% CALCULATIONS:
% -------------
% -------------------------------------------------------------
% PART 1: Redo Project 2 with REAL products after the combustor
% --------------------------------------------------------------
% Assume: complete combustion of Jet A


% PROJECT 2, PART 2: Spool Speed vs. T0,P0,Ma,U,mdot
% --------------------------------------------------
% NOTE:
% (i) UNITS: all pressures in KPa, temperatures in C, areas in in^2
% (ii) subscript m = measured temp before Recovery factor adjsutment
% (ii) no subscript = actual stagnation 
% (ii) subscript 0 = actual stagnaton 

% Station 2 (still need to get gamma at T2m since can't assume constant specific heat)
P02 = linspace(Patm, Patm, length(rpm))';               % P02 = P01 [ASSUME: isentropic from 0->2]
P2 = P02 - dP2;                                         % P2 [absolute] 
[~,~,k2] = sp_heats(T2m);

M2 = sym('m',[length(rpm),1]);
assume(M2, 'real');
assumeAlso(M2 > 0); %Select the subsonic solution. Set M > 1 if want supersonic solution.
Ma2 = vpasolve(P02./P2 ==(1+ M2.^2.*((k2-1)./2)).^(k2./(k2-1)), M2); % Solve for Ma2 [ASSUME: isentropic from 0->2]
Ma2 = struct(Ma2);
Ma2 = struct2array(Ma2)'; %Ma is returned as a symbolic variable which doesn't work in polyval in sp_heats
Ma2 = abs(Ma2);           %This seems kind of sketchy, but looking at the results seems reasonable. Tested positive solution and it is valid, given M2 is squared in all terms.
Ma2 = double(Ma2);
T2  = recoveryFactor(T2m,Ma2,RFcross);
T02 = T2.*(1 + (Ma2.^2).*((k2-1)./2));
U2 = Ma2.*sqrt(k2.*R.*T2);
MFP2 = findMFP(Ma2,k2);
             
% Mass Flow of Air and Fuel & Air/Fuel Ratio Plots
mdot_air = ( MFP2.*A2.*P02 )./sqrt(R.*T02);                 % USE: MFP = (mdot/A)*(sqrt(R*T0)/P0)
mdot_total = mdot_air + mdot_fuel;
AF = mdot_air ./ mdot_fuel;                                 % air fuel ratio
phi = AFs./AF;

% Station 3 (compressor exit)
[T3,T03,Ma3,U3,k3] = findMaTemps(mdot_total,A3,T3m,P03,RFcross);

% CODE CHANGES FROM PROJECT # 2 HERE
% ------------------------------------
% Station 4 (after combustor, before bend to turbine inlet)
P04_guess = P4;                                         % P04 = P4 [ASSUME: Ma~<0.1 and delta(P0)~5%, therefore P04~P4] see LEC 6 page 23
[T4,T04,Ma4,U4,k4] = findMaTemps_mix(mdot_total,A4,T4m,P04_guess,RFcross,AF);
P04 = P4.*(1+ Ma4.^2.*((k4-1)./2)).^(k4./(k4-1));       


% Station 5 (turbine outlet)                                  
[T5,T05,Ma5,U5,k5] = findMaTemps_mix(mdot_total,A5,T5m,P05,RFaxial,AF);

% Station 8 (nozzle exit)            
[T8,T08,Ma8,U8,k8] = findMaTemps_mix(mdot_total,A8,T8m,P08,RFcross,AF);
P8 = P08./(1+ Ma8.^2.*((k8-1)./2)).^(k8./(k8-1));       

% Station 1 (can't be found first)
T01 = T02;                                              % [ASSUME: Isentropic from 0-->2, therefore also isentropic 1-->2]
P01 = P02;                                              % [ASSUME: small pressure change in combustor.]

% Find Net Thrust
% INLET: state 1 (beg. of nozzle), OUTLET: state 8 (outlet of nozzle)
U1 = findU(T01,P01,mdot_air,A1);
u_in = U1;
u_out = U8;
Fthrust = (mdot_total .* u_out) - (mdot_air .* u_in);
ST = Fthrust./mdot_total;
Qdot = mdot_fuel*LHV(1:2);                  
TSFC = mdot_fuel./Fthrust;

% PART 2: PLOTS
% -------------
% Stagnation Temperature Plots
figure(1)
%subplot(2,2,1)
plot(rpm./1000,T02,'b',rpm./1000,T03,'c',rpm./1000,T04,'m',rpm./1000,T05,'r',rpm./1000,T08,'g');
title('Part 2: Station Stagnation Temperature vs Spool Speed');
ylabel('Stagnation Temperature [K]');
xlabel('Spool Speed [rpm x1000]');
legend('T02','T03','T04','T05','T08', 'Location', 'East');
if plotfixing; plotfixer; end

% Stagnation Pressure Plots
figure(2)
%subplot(2,2,2)
plot(rpm./1000,P02./KPA_TO_PA,'b',rpm./1000,P03./KPA_TO_PA,'c',rpm./1000,P04./KPA_TO_PA,'m',rpm./1000,P05./KPA_TO_PA,'r',rpm./1000,P08./KPA_TO_PA,'g');
title('Part 2: Station Stagnation Pressure vs Spool Speed');
ylabel('Stagnation Pressure [kPa]'); 
xlabel('Spool Speed [rpm x1000]');
legend('P02','P03','P4','P05','P08', 'Location', 'East');
if plotfixing; plotfixer; end


% Mach Plots
figure(3)
%subplot(2,2,3)
plot(rpm./1000,Ma2,'b',rpm./1000,Ma3,'c',rpm./1000,Ma4,'m',rpm./1000,Ma5,'r',rpm./1000,Ma8,'g');
title('Part 2: Station Mach Number vs Spool Speed');
ylabel('Mach number');
xlabel('Spool Speed [rpm x1000]');
legend('Ma2','Ma3','Ma4','Ma5','Ma8', 'Location', 'NorthWest');
if plotfixing; plotfixer; end


% Velocity Plots
figure(4)
%subplot(2,2,4)
plot(rpm./1000,U2,'b',rpm./1000,U3,'c',rpm./1000,U4,'m',rpm./1000,U5,'r',rpm./1000,U8,'g');
title('Part 2: Station Velocity vs Spool Speed');
ylabel('Velocity [m/s]');
xlabel('Spool Speed [rpm x1000]');
legend('U02','U03','U04','U05','U08','Location', 'NorthWest');
if plotfixing; plotfixer; end


% Mass Flow Plots
figure(5)
plot(rpm./1000,mdot_air.*KG_TO_G,'b', rpm./1000,mdot_fuel.*KG_TO_G,'c',rpm./1000,AF,'m');
title('Part 2: Mass flow Rates vs Spool Speed');
ylabel('Mass FLow Rate [g/s]');
xlabel('Spool Speed [rpm x1000]');
legend('Mass Flow Air','Mass Flow Fuel','Air-Fuel Ratio', 'Location', 'NorthWest');
if plotfixing; plotfixer; end


% Thrust Flow Plots
figure(6)
plot(rpm./1000,Fthrust,'b', rpm./1000,Fthrust_meas,'c');
title('Part 2: Spool Speed vs Net Thrust and Measured Thrust ');
ylabel('Thrust [N]');
xlabel('Spool Speed [rpm x1000]');
legend('Fthrust (calculated)','Fthrust (measured)', 'Location', 'NorthWest');
if plotfixing; plotfixer; end


% PROJECT 2 PART 3: Spool Speed vs. ST TSFC, Thermal Efficiency
% -------------------------------------------------------------
Wnet = 0.5.*( mdot_total.*u_out.^2 - mdot_air.*u_in.^2 );
eta_thermal = zeros(size(Qdot));
for i = 1:size(Qdot, 2);
    eta_thermal(:,i) = Wnet./Qdot(:,i);
end

% Thrust Flow Plots
figure(7)
%subplot(2,2,1)
plot(rpm./1000,ST,'b')
title('Part 3: Specific Thrust vs Spool Speed ');
ylabel('Thrust [N]');
xlabel('Spool Speed [rpm x1000]');
if plotfixing; plotfixer; end


figure(8)
%subplot(2,2,2)
plot(rpm./1000,TSFC.*1000,'m')
title('Part 3: Thrust Specific Fuel Consumption vs Spool Speed');
ylabel('Thrust [g*N/s]');
xlabel('Spool Speed [rpm x1000]');
if plotfixing; plotfixer; end


figure(9)
%subplot(2,2,3)
plot(rpm./1000,eta_thermal);
title('Part 3: Thermal Efficiency vs Spool Speed ');
ylabel('Thermal Effeciency');
xlabel('Spool Speed [rpm x1000]');
legend('Jet-A','Dodecane');
if plotfixing; plotfixer; end


% PROJECT 2 PART 4: Spool Speed vs. Power Compressor, Power Turbine 
% --------------------------------------------------------
% Adiabatic Efficiencies
% COMPRESSOR (2->3): total-to-total, Win,isentropic / Win,actual
% COMBUSTOR (3->4)
% TURBINE (4->5): total-to-total, assume exhaust KE is not a loss, aeropropulsive
% NOZZLE (5->8): total-to-static

Pin_compressor = find_dh_mix( T02,T03,phi );
Pout_turbine_project3 = find_dh_mix( T04,T05,phi );

% Compressor Outlet
T03s = applyIsentropicTempVar( T02, P03./P02 );       % tt

% Turbine
T05s = applyIsentropicTempVar_mix( T04, P05./P04, phi );   % tt

% Nozzle
T8s = applyIsentropicTempVar_mix( T05, P8./P05, phi );     % ts

f_temp = @(x) sp_heats(x); %Note: Cp is the first element in the return of sp_heats and this returns just cp as a 1x1 double
f_temp_mix = @(x) sp_heats_mix(x,phi);
n_compressor = zeros(length(T02),1);
n_turbine =    zeros(length(T02),1);
n_nozzle =     zeros(length(T02),1);
for i = 1:length(T02)    
    n_compressor(i) = integral( f_temp,T02(i),T03s(i) )/integral( f_temp,T02(i),T03(i) );
    n_turbine(i) =    integral( f_temp_mix,T05(i),T04(i) )/integral( f_temp_mix,T05s(i),T04(i) );
    n_nozzle(i) = (.5.*U8(i).^2)./integral( f_temp_mix,T8s(i),T05(i) );
end

P0ratio_combustor = P04./P03;

% ------------------------------------------
% PART 1: Turbine Power Mixed vs. Project #2 
% ------------------------------------------
Pout_turbine_project2 = me140_project2_Pturbine();
Pout_turbine_300K = me140_project2_300K();

% ---------------------------------
% PART 2: hf & Adiabatic Flame Temp
% ---------------------------------
TR = 300; % [K]  
phi_range = linspace(0.05, 0.65)';
[ hf_jetA, TP_part2 ] = findAdiabaticFlameTemp( hf_co2, hf_h2o, hf_n2, hf_o2, LHV(1:2), Mfuel(1), TR , phi_range );

% -----------------------------------------------------
% PART 3: Turbine Power using Adiabatic Burned Gas Temp
% -----------------------------------------------------
% Find Adiabatic Flame Temp, TP = T04 (using T03)
[~, TP_part3] = findAdiabaticFlameTemp( hf_co2, hf_h2o, hf_n2, hf_o2, LHV(1:2), Mfuel(1), T03, phi );

% Find T05s (using T04 = Tp)
T05s_new = applyIsentropicTempVar( TP_part3, P05./P04 ); 

% Find Isentropic Turbine Work: Wturbine,isentropic = -mdot*Cp*(T05s-T04);
Pout_turbine_project3 = find_dh_mix( TP_part3,T05s_new,phi );

% ----------------------------------------
% PART 4: ENGINE PERFORMANCE IMPROVEMENTS
% ----------------------------------------
% (1) EFFECT OF INCREASING LOWER HEATING VALUE (LHV)
% Pout_turbine_project3 = find_dh_mix( T04,T05,phi );

% PART 1 PLOTS
% ------------
% Power Plots
figure(10)
plot(rpm ./1000, Pin_compressor./1000, 'r', rpm./1000, abs(Pout_turbine_project3)./1000, 'b');
xlabel('Spool Speed [rpm x1000]');
ylabel('Power [kW]');
title('Power into Compressor and Power Out of Turbine vs. Spool Spped');
legend('Power into Compressor', 'Power Out of Turbine', 'Location', 'NorthWest');
if plotfixing; plotfixer; end

% Adiabatic Efficiency Plots
figure(11)
plot(rpm./1000, n_compressor, 'r', rpm./1000, n_turbine, 'b', rpm./1000, n_nozzle, 'm');
title('Efficiency of Compressor, Turbine, Nozzle vs Spool Speed');
xlabel('Spool Speed [rpm x1000]');
ylabel('Efficiency');
legend('Compressor Efficiency', 'Turbine Efficiency', 'Nozzle Efficiency');
if plotfixing; plotfixer; end

% Stagnation Pressure Ratio Across Combustor
figure(12)
plot(rpm./1000, P0ratio_combustor, 'g');
title('Stagnation Pressure Across Combustor vs Spool Speed');
xlabel('Spool Speed [rpm x1000]');
ylabel('Stagnation Pressure Ratio, P04/P03');
if plotfixing; plotfixer; end

% Turbine Power, 3 Ways
figure(13)
plot(rpm./1000, abs(Pout_turbine_project3)./1000, 'g',rpm./1000,abs(Pout_turbine_project2)./1000,'c',rpm./1000,abs(Pout_turbine_300K)./1000,'m');
title('Turbine Power vs Spool Speed');
xlabel('Spool Speed [rpm x1000]');
ylabel('Turbine Power [kW]');
legend('Project 3: accounted for air/fuel mixing','Project 2: Variable Cp','Project 2: Assuming constant Cp = Cp(300K)');
if plotfixing; plotfixer; end

% PART 2 PLOTS
% ------------
% Adiabatic Flame Temp as a Function of Phi
figure(14)
plot(phi_range, TP_part2,'b');
xlabel('Equivalence Ratio, \phi');
ylabel('Adiabatic Flame Temperature, Tp [K]');
title('Adiabatic Flame Temp as a Function of Phi');
if plotfixing; plotfixer; end
