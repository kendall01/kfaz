function [ Pout_turbine ] = me140_project2_300K(~)
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
plotfixing = 0;       % Boolean for whether or not to run plot fixer

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
dP2 = me140_project2_data(8)' .* KPA_TO_PA;            % dP2 = Pstag - Pstatic 
P03 = me140_project2_data(9)' .* KPA_TO_PA + Patm ;    % P03 [absolute] [Pa]
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
latent_heat_vap= [42800*KJ_TO_J 44560*KJ_TO_J];        % latent heat of vaporization
LHV = [42800*KJ_TO_J, 44560*KJ_TO_J];                  % vapor, J/kg
hf0 = [0 -1712];                                       % vapor, enthalpy of formation J/kg (replace 0 later)
AFs = 14.43;                    % stoichiometric Air-Fuel-Ratio, found from balancing equation (LEC 5, slides 3-4)


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
[~,~,k2] = sp_heats_300K(T2m);

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
Qdot = mdot_fuel.*LHV(1);                  
TSFC = mdot_fuel./Fthrust;

% PROJECT 2 PART 3: Spool Speed vs. ST TSFC, Thermal Efficiency
% -------------------------------------------------------------
Wnet = 0.5.*( mdot_total.*u_out.^2 - mdot_air.*u_in.^2 );
eta_thermal = Wnet./Qdot;


% PROJECT 2 PART 4: Spool Speed vs. Power Compressor, Power Turbine 
% --------------------------------------------------------
% Adiabatic Efficiencies
% COMPRESSOR (2->3): total-to-total, Win,isentropic / Win,actual
% TURBINE (4->5): total-to-total, assume exhaust KE is not a loss, aeropropulsive
% NOZZLE (5->8): total-to-static

Pin_compressor = find_dh_mix( T02,T03,AFs./AF );
Pout_turbine = find_dh_mix( T04,T05,AFs./AF );
end

