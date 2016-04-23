% Project4.m
% 4-22-16 - Created Jon Renslo
% Script for Project 4 in ME140. Fuel Cells. 
close all; clear; clc;
%% Part 1
% find first law eta_max for a PEM fuel cell 
npts = 100;
T = linspace(25+273,1273,npts);
lambda = 2; % 100% excess air
HHV_h2 = 141800; %kj/kg
LHV_h2 = 120000; %kj/kg
MMh2 = 2.02; %g/mol
gToKg = 1000;
P = 101.3e3; %Pa, Preact = Pprod = Patm

% assume 1 mol of h2 combusted --> 4.76/2 mol air
% 139 g total
molh2 = 1;
massh2 = molh2*MMh2/gToKg;
molair = 4.76*lambda/2;
molo2_rxn = molair/4.76;
moln2 = molair*3.76/4.76;

molh2o = molh2;
molo2_prod = 0.5*(lambda-molh2)*molo2_rxn;

% check mass balance
% mreact = massh2 + molair*28.85/1000
% mprod = molo2_prod*32/1000 + moln2*28/1000 + molh2o*18.02/1000

% gives in Joules, as we multiply by moles
gprod_LHV = gEng(T,P,'h2ovap',molh2o) + gEng(T,P,'o2',molo2_prod) + gEng(T,P,'n2',moln2);
gprod_HHV = gEng(T,P,'h2o',molh2o) + gEng(T,P,'o2',molo2_prod) + gEng(T,P,'n2',moln2);
greact = gEng(T,P,'h2',molh2) + gEng(T,P,'o2',molo2_rxn) + gEng(T,P,'n2',moln2);

delG_hhv = gprod_HHV-greact;
delG_lhv = gprod_LHV-greact;
eta_hhv = -delG_hhv/(HHV_h2*massh2*1000);
eta_lhv = -delG_lhv/(LHV_h2*massh2*1000);
%todo find g actual based on liquid water mixture
figure();
plot(T,eta_hhv,T,eta_lhv);
legend('\eta_{HHV}','\eta_{LHV}');
xlabel('Temperature K');
ylabel('Maximum 1st Law Efficiency');
plotfixer();


