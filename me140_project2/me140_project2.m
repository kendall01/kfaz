% ME 140 PROJECT 2: SR30 TURBOJET ENGINE
% FRANKIE WILLCOX, KENDALL FAGAN, ZACH NEVILLS, ANTOINE SCREVE
% -------------------------------------------------------------

clc;
clear all;
clearvars;

% UNIT CONVERSIONS
% ----------------
C_TO_K = 273.15;
LBF_TO_N = 4.44822;
IN2_TO_M2 = 6.4516*10^-4;
KPA_TO_PA = 10^3;
KG_TO_G = 10^3;
R = 287.058;                                           % J/(kg*K)

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
dP2 = me140_project2_data(8)' .* KPA_TO_PA;    % dP2 = Pstag - Pstatic 
P03 = me140_project2_data(9)' .* KPA_TO_PA + Patm ;    % P03 [absolute] [Pa]
P4 =  me140_project2_data(10)' .* KPA_TO_PA + Patm;    % P4 (static) [absolute] 
P05 = me140_project2_data(11)' .* KPA_TO_PA + Patm;    % P05 [absolute]   
P08 = me140_project2_data(12)' .* KPA_TO_PA + Patm;    % P08 [absolute]   


% GIVEN:
% ------
RFcross = linspace(.68, .68, length(rpm))';            % Reference Factor cross-section
RFaxial = linspace(.86, .86, length(rpm))';            % Reference Factor axial

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
% PART 2: Spool Speed vs. T0,P0,Ma,U,mdot
% -------------------------------------------------------------
% NOTE:
% (i) UNITS: all pressures in KPa, temperatures in C, areas in in^2
% (ii) subscript m = measured temp before Recovery factor adjsutment
% (ii) no subscript = actual stagnation 
% (ii) subscript 0 = actual stagnaton 

% Station 2 (still need to get gamma at T2m since can't assume constant specific heat)
P02 = linspace(Patm, Patm, length(rpm))';               % P02 = P01 [ASSUME: isentropic from 0->2]
P2 = P02 - dP2;                                         % P2 [absolute] 
k2 = sp_heats(T2m);

M2 = sym('m',[length(rpm),1]);
assume(M2, 'real');
assumeAlso(M2 > 0); %Select the subsonic solution. Set M > 1 if want supersonic solution.
Ma2 = vpasolve(P02./P2 ==(1+ M2.^2.*((k2-1)./2)).^(k2./(k2-1)), M2); % Solve for Ma2 [ASSUME: isentropic from 0->2]
Ma2 = struct(Ma2);
Ma2 = struct2array(Ma2)'; %Ma is returned as a symbolic variable which doesn't work in polyval in sp_heats
Ma2 = abs(Ma2);     %This seems kind of sketchy, but looking at the results seems reasonable. Tested positive solution and it is valid, given M2 is squared in all terms.
Ma2 = double(Ma2);
T2  = recoveryFactor(T2m,Ma2,RFcross);
T02 = T2.*(1 + (Ma2.^2).*((k2-1)./2));
U2 = Ma2.*sqrt(k2.*R.*T2);
MFP2 = findMFP(Ma2,k2);
             
% Mass Flow of Air and Fuel & Air/Fuel Ratio Plots
mdot_air = ( MFP2.*A2.*P02 )./sqrt(R.*T02);                 % USE: MFP = (mdot/A)*(sqrt(R*T0)/P0)
mdot_total = mdot_air + mdot_fuel;

% Station 3 (compressor exit)
[T3,T03,Ma3,U3,k3] = findMaTemps(mdot_total,A3,T3m,P03,RFcross);


% Station 4 (before bend to turbine inlet)
P04_guess = P4;                                 % P04 = P4 [ASSUME: Ma~<0.1 and delta(P0)~5%, therefore P04~P4] see LEC 6 page 23
[T4,T04,Ma4,U4,k4] = findMaTemps(mdot_total,A4,T4m,P04_guess,RFcross);
P04 = P4.*(1+ Ma4.^2.*((k4-1)./2)).^(k4./(k4-1));     % Use P4,T4,T04 to find REAL P04

% Station 5 (turbine outlet)                                  
[T5,T05,Ma5,U5,k5] = findMaTemps(mdot_total,A5,T5m,P05,RFaxial);

% Station 8 (nozzle exit)            
[T8,T08,Ma8,U8,k8] = findMaTemps(mdot_total,A8,T8m,P08,RFcross);

% Station 1 (can't be found first)
% Process:
% - find Win turbine using 4->5 ??
% - entropy balance with fuel in, work in, mass in, mass out ??


% PART 2: PLOTS
% -------------
% Stagnation Temperature Plots
figure(1)
plot(rpm,T02,'b',rpm,T03,'c',rpm,T04,'m',rpm,T05,'r',rpm,T08,'g');
title('Spool Speed vs Station Stagnation Temperature');
ylabel('Stagnation Temperature [K]');
xlabel('Spool Speed [rpm]');
legend('T02','T03','T04','T05','T08');

% Stagnation Pressure Plots
figure(2)
plot(rpm,P02,'b',rpm,P03,'c',rpm,P04,'m',rpm,P05,'r',rpm,P08,'g');
title('Spool Speed vs Station Stagnation Pressure');
ylabel('Stagnation Pressure [Pa]'); 
xlabel('Spool Speed [rpm]');
legend('P02','P03','P4','P05','P08');

% Mach Plots
figure(3);
plot(rpm,Ma2,'b',rpm,Ma3,'c',rpm,Ma4,'m',rpm,Ma5,'r',rpm,Ma8,'g');
title('Spool Speed vs Station Mach Number');
ylabel('Mach number');
xlabel('Spool Speed [rpm]');
legend('Ma2','Ma3','Ma4','Ma5','Ma8');

% Velocity Plots
figure(4);
plot(rpm,U2,'b',rpm,U3,'c',rpm,U4,'m',rpm,U5,'r',rpm,U8,'g');
title('Spool Speed vs Station Velocity');
ylabel('Velocity [m/s]');
xlabel('Spool Speed [rpm]');
legend('U02','U03','U04','U05','U08');

% % Mass Flow Plots
figure(5)
plot(rpm,mdot_air.*KG_TO_G,'b', rpm,mdot_fuel.*KG_TO_G,'c',rpm,air_fuel_ratio,'m');
title('Spool Speed vs Mass flow Rates');
ylabel('Mass FLow Rate [g/s]');
xlabel('Spool Speed [rpm]');
legend('Mass Flow Air','Mass Flow Fuel','Air-Fuel Ratio');

% -------------------------------------------------------------
% PART 3: Spool Speed vs. ST TSFC, Thermal Efficiency
% -------------------------------------------------------------
% eta_thermal = 

% -------------------------------------------------------------
% PART 4: Spool Speed vs. Power Compressor, Power Turbine and OTHER CALCS
% -------------------------------------------------------------
