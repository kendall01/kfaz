% ME 140 PROJECT 2: SR30 TURBOJET ENGINE
% FRANKIE WILLCOX, KENDALL FAGAN, ZACH NEVILLS, ANTOINE SCREVE
% -------------------------------------------------------------

clc;
clear all;
clearvars;
% DATA:
% ------
rpm = me140_project2_data(1)';
mdot_fuel = me140_project2_data(13).*1000;              % mdot_fuel,[g/s]
Fthrust_meas = me140_project2_data(14).*4.44822;        % Fthrust, [N]
Toil = me140_project2_data(7);

% GIVEN:
% ------
rho2 = 0.074887*(1/12)^3; % [lb/in^3]
RFcross = linspace(.68, .68, length(rpm))'; %Reference Factor cross-section
RFaxial = linspace(.86, .86, length(rpm))'; %Reference Factor axial
R = 287; 
Patm = 101.3;

% Areas [in^2]
A1 = 27.3; % Flow area at inlet of bellmouth
A2 = 6.4;  % Effective flow area at pitot-static probe
A3 = 9.0;  % Area of compressor exit
A4 = 7.2;  % Area before bend to turbine inlet
A5 = 4.7;  % Area of turbine outlet
A8 = 3.87; % Area of nozzle exit

% CALCULATIONS:
% -------------


% -------------------------------------------------------------
% PART 1: Tables of Data
% -------------------------------------------------------------

% All tables uploaded to google doc and made using excel

% -------------------------------------------------------------
% PART 2: Spool Speed vs. T0,P0,Ma,U,mdot
% -------------------------------------------------------------
% NOTE:
% (i) UNITS: all pressures in KPa, temperatures in C, areas in in^2
% (ii) subscript m = measured temp before Recovery factor adjsutment
% (ii) no subscript = actual stagnation 
% (ii) subscript 0 = actual stagnaton 

% TO-DO's:
% ---------
% - check where we assumed P0=P
% - check that we want to use Cp(T) instead of Cp(T0) in findMaTemps
% function, but ok to initialize as Cp(Tm)

% Station 2 (still need to get gamma at T2m since can't assume constant specific heat)
dP2 = me140_project2_data(8)';                   % dP2 = Pstag - Pstatic 
P02 = linspace(Patm, Patm, length(rpm))';       % P02 = P01 [ASSUME: isentropic from 0->2]
P2 = P02 - dP2;                                 % P2 [absolute] 
T2m = me140_project2_data(2)';                   % T2m, cross-flow
k2 = sp_heats(T2m);

M2 = sym('m',[length(rpm),1]);
assume(M2, 'real');
assumeAlso(M2 > 0); %Select the subsonic solution. Set M > 1 if want supersonic solution.
Ma2 = vpasolve(P02./P2 ==(1+ M2.^2.*((k2-1)./2)).^(k2./(k2-1)), M2); % Solve for Ma2 [ASSUME: isentropic from 0->2]
Ma2 = struct(Ma2);
Ma2 = struct2array(Ma2)'; %Ma is returned as a symbolic variable which doesn't work in polyval in sp_heats
Ma2 = abs(Ma2);     %This seems kind of sketchy, but looking at the results seems reasonable. Tested positive solution and it is valid, given M2 is squared in all terms.
T2  = recoveryFactor(T2m,Ma2,RFcross);
T02 = T2.*(1 + (Ma2.^2).*((k2-1)./2));
U2 = Ma2.*sqrt(k2.*R.*T2);
MFP2 = findMFP(Ma2,k2);
             
% Mass Flow of Air and Fuel & Air/Fuel Ratio Plots
mdot_air = ( MFP2*A2*P02 )/sqrt(R*T02);                 % USE: MFP = (mdot/A)*(sqrt(R*T0)/P0)
mdot_total = mdot_air + mdot_fuel;

% Station 3 (compressor exit)
P03 = me140_project2_data(9)+101.3;             % P03 [absolute] 
T3m = me140_project2_data(3);                   % T3m, cross-flow 
[T3,T03,Ma3,U3,k3] = findMaTemps(mdot_total,A3,P03,T3m,RFcross);

% Station 4 (before bend to turbine inlet)
P4 = me140_project2_data(10)+101.3;             % P4 (static) [absolute] 
P04_guess = P4;                                 % P04 = P4 [ASSUME: Ma~<0.1 and delta(P0)~5%, therefore P04~P4] see LEC 6 page 23
T4m = me140_project2_data(4);                   % T4m, cross-flow
[T4,T04,Ma4,U4,k4] = findMaTEmps(mdot_total,P04_guess,T4m,RFcross);
P04 = P4*(1+ Ma4^2*((k4-1)/2))^(k4/(k4-1));     % Use P4,T4,T04 to find REAL P04

% Station 5 (turbine outlet)
P05 = me140_project2_data(11)+101.3;            % P05 [absolute]                                    
T5m = me140_project2_data(5);                   % T5m, axial-flow 
[T5,T05,Ma5,U5,k5] = findMaTemps(mdot_total,P05,T5m,RFaxial);

% Station 8 (nozzle exit)
P08 = me140_project2_data(12)+101.3;            % P08 [absolute]               
T8m = me140_project2_data(6);                   % T8m, cross-flow
[T8,T08,Ma8,U8,k8] = findMaTemps(mdot_total,P08,T8m,RFcross);

% Station 1 (can't be found first)
% Process:
% - find Win turbine using 4->5 ??
% - entropy balance with fuel in, work in, mass in, mass out ??


% PART 2: PLOTS
% -------------
% Stagnation Temperature Plots
figure(1)
plot(T02,rpm,'b',T03,rpm,'c',T04,rpm,'m',T05,rpm,'r',T08,rpm,'g');
title('Spool Speed vs Station Stagnation Temperature');
xlabel('Stagnation Temperature (Celsius)');
ylabel('Spool Speed (rpm)');
legend('T02','T03','T04','T05','T08');

% Stagnation Pressure Plots
figure(2)
plot(P02,rpm,'b',P03,rpm,'c',P04,rpm,'m',P05,rpm,'r',P08,rpm,'g');
title('Spool Speed vs Station Stagnation Pressure');
xlabel('Stagnation Pressure (kPa)'); 
ylabel('Spool Speed (rpm)');
legend('P02','P03','P4','P05','P08');

% Mach Plots
figure(3);
plot(Ma2,rpm,'b',Ma3,rpm,'c',Ma4,rpm,'m',Ma5,rpm,'r',Ma8,rpm,'g');
title('Spool Speed vs Station Mach Number');
xlabel('Mach number');
ylabel('Spool Speed (rpm)');
legend('Ma2','Ma3','Ma4','Ma5','Ma8');

% Velocity Plots
figure(4);
plot(U2,rpm,'b',U3,rpm,'c',U4,rpm,'m',U5,rpm,'r',U8,rpm,'g');
title('Spool Speed vs Station Velocity');
xlabel('Velocity (m/s)');
ylabel('Spool Speed (rpm)');
legend('U02','U03','U04','U05','U08');

% Mass Flow Plots
figure(5)
plot(mdot_air,rpm,'b', mdot_fuel,rpm,'c',air_fuel_ratio,rpm,'m');
title('Spool Speed vs Mass flow Rates');
xlabel('Velocity (m/s)');
ylabel('Spool Speed (rpm)');
legend('Mass Flow Air','Mass flow Fuel','Air-Fuel Ratio');

% -------------------------------------------------------------
% PART 3: Spool Speed vs. ST TSFC, Thermal Efficiency
% -------------------------------------------------------------
% eta_thermal = 

% -------------------------------------------------------------
% PART 4: Spool Speed vs. Power Compressor, Power Turbine and OTHER CALCS
% -------------------------------------------------------------
